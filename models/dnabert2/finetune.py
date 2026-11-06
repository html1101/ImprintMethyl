# DNABERT-2 first-pass fine-tune: sequence-only (A) vs sequence + parent-of-origin (B).
# Input: 1 kb mm9 reference window centred on a CpG. Target: per-allele methylation fraction.
import os, sys, io, csv, gzip, json, math, random, tarfile, time, zipfile
sys.modules["triton"] = None  # DNABERT-2's Triton attention does not compile here; use plain attention
import torch, torch.nn as nn
from transformers import AutoTokenizer, AutoModel

BASE = os.environ.get("IMPRINT_BASE", "/ocean/projects/cis260266p/rpandit1/imprintmethyl")
OUT = os.environ.get("OUT", BASE + "/runs/run_" + os.environ.get("SLURM_JOB_ID", "local"))
os.makedirs(OUT, exist_ok=True)
WIN = int(os.environ.get("WIN", "1000"))
TRAIN_PCT = int(os.environ.get("TRAIN_PCT", "5"))
MAXTRAIN = int(os.environ.get("MAXTRAIN", "400000"))
MAXTEST = int(os.environ.get("MAXTEST", "40000"))
MINCOV = int(os.environ.get("MINCOV", "5"))
EPOCHS = float(os.environ.get("EPOCHS", "1"))
BATCH = int(os.environ.get("BATCH", "64"))
TEST_CHROM = "chr19"
AUTOS = {"chr%d" % i for i in range(1, 20)}
random.seed(0); torch.manual_seed(0)

def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)

# ---- imprinted DMRs (Xie 2012 Table S2, mm9) ----
dmr_by_chr = {}
for r in csv.DictReader(open(BASE + "/ref/xie2012_imprinted_dmrs_mm9.csv")):
    dmr_by_chr.setdefault(r["chrom"], []).append((int(r["start_mm9"]), int(r["end_mm9"]), r["methylated_parent"], r["locus"]))
def in_dmr(c, p):
    for s, e, par, loc in dmr_by_chr.get(c, ()):
        if s <= p <= e:
            return (par, loc)
    return None

# ---- stream Sarah's table, merge strands per CpG (watson pos = C; crick pos - 1 = same CpG) ----
log("reading methylation table")
zf = zipfile.ZipFile(BASE + "/processed/processed.zip")
fh = io.TextIOWrapper(gzip.GzipFile(fileobj=zf.open("processed/methylation_cpg_counts_mm9.tsv.gz")))
next(fh)
agg = {}
n = 0
for line in fh:
    n += 1
    f = line.split("\t")
    c = f[2]
    if c not in AUTOS:
        continue
    pos = int(f[3])
    cpg = pos if f[1] == "watson" else pos - 1
    d = in_dmr(c, cpg)
    keep = (c == TEST_CHROM) or (d is not None) or (((cpg * 2654435761) & 0xffffffff) % 100 < TRAIN_PCT)
    if not keep:
        continue
    k = (f[0], c, cpg)
    v = agg.get(k)
    if v is None:
        v = agg[k] = [0, 0, 0, 0]
    v[0] += int(f[13]); v[1] += int(f[14]); v[2] += int(f[17]); v[3] += int(f[18])
log("rows read", n, "CpG x cross kept", len(agg))

# ---- examples: two per CpG x cross (maternal, paternal) ----
cpgs = {}          # (chrom, cpg) -> index into window list
rec = []           # (cpg_idx, cross, parent 0=mat 1=pat, frac, split, dmr_parent, dmr_locus)
for (cross, c, cpg), (mm, mu, pm, pu) in agg.items():
    if mm + mu < MINCOV or pm + pu < MINCOV:
        continue
    d = in_dmr(c, cpg)
    split = "dmr" if d else ("test" if c == TEST_CHROM else "train")
    idx = cpgs.setdefault((c, cpg), len(cpgs))
    rec.append((idx, cross, 0, mm / (mm + mu), split, d[0] if d else "", d[1] if d else ""))
    rec.append((idx, cross, 1, pm / (pm + pu), split, d[0] if d else "", d[1] if d else ""))
del agg
log("examples", len(rec), "unique CpGs", len(cpgs))

# ---- sequence windows from mm9 ----
log("extracting windows")
windows = [None] * len(cpgs)
by_chr = {}
for (c, cpg), i in cpgs.items():
    by_chr.setdefault(c, []).append((cpg, i))
tf = tarfile.open(BASE + "/ref/chromFa.tar.gz")
names = {m.name.split("/")[-1]: m for m in tf.getmembers()}
cg_ok = cg_n = 0
for c, lst in by_chr.items():
    raw = tf.extractfile(names[c + ".fa"]).read().decode()
    seq = "".join(raw.split("\n")[1:]).upper()
    for cpg, i in lst:
        c0 = cpg - 1
        cg_n += 1; cg_ok += seq[c0:c0 + 2] == "CG"
        s = seq[c0 - WIN // 2: c0 + WIN // 2]
        if len(s) == WIN and s.count("N") < WIN // 10:
            windows[i] = s
    del seq, raw
log("CpG check: %d/%d windows centred on CG (%.2f%%)" % (cg_ok, cg_n, 100.0 * cg_ok / max(cg_n, 1)))
rec = [r for r in rec if windows[r[0]] is not None]
train = [r for r in rec if r[4] == "train"]; random.shuffle(train); train = train[:MAXTRAIN]
test = [r for r in rec if r[4] == "test"]; random.shuffle(test); test = test[:MAXTEST]
dmr = [r for r in rec if r[4] == "dmr"]
log("train %d  test(chr19) %d  dmr %d" % (len(train), len(test), len(dmr)))

# ---- tokenize unique windows once ----
name = "zhihan1996/DNABERT-2-117M"
tok = AutoTokenizer.from_pretrained(name, trust_remote_code=True)
used = sorted({r[0] for r in train + test + dmr})
pos_of = {j: k for k, j in enumerate(used)}
enc = tok([windows[j] for j in used], padding="longest", truncation=True, max_length=512, return_tensors="pt")
IDS, AM = enc["input_ids"], enc["attention_mask"]
log("tokenized", IDS.shape)

def batches(rs, bs, shuffle):
    rs = list(rs)
    if shuffle:
        random.shuffle(rs)
    for i in range(0, len(rs), bs):
        b = rs[i:i + bs]
        k = torch.tensor([pos_of[r[0]] for r in b])
        yield (IDS[k].cuda(non_blocking=True), AM[k].cuda(non_blocking=True),
               torch.tensor([r[2] for r in b]).cuda(), torch.tensor([r[3] for r in b], dtype=torch.float32).cuda(), b)

class Net(nn.Module):
    def __init__(self, use_parent):
        super().__init__()
        self.b = AutoModel.from_pretrained(name, trust_remote_code=True)
        self.use = use_parent
        self.pe = nn.Embedding(2, 32)
        self.head = nn.Sequential(nn.Linear(768 + (32 if use_parent else 0), 256), nn.GELU(), nn.Linear(256, 1))
    def forward(self, ids, am, parent):
        h = self.b(ids, attention_mask=am)[0]
        x = ((h * am.unsqueeze(-1)).sum(1) / am.sum(1, keepdim=True)).float()
        if self.use:
            x = torch.cat([x, self.pe(parent)], -1)
        return torch.sigmoid(self.head(x)).squeeze(-1)

def pearson(a, b):
    a = torch.tensor(a); b = torch.tensor(b)
    a = a - a.mean(); b = b - b.mean()
    return float((a * b).sum() / (a.norm() * b.norm() + 1e-12))

def evaluate(model, rs):
    model.eval(); preds = []
    with torch.no_grad():
        for ids, am, par, y, b in batches(rs, 256, False):
            with torch.autocast("cuda", dtype=torch.bfloat16):
                p = model(ids, am, par)
            preds += p.float().cpu().tolist()
    model.train()
    ys = [r[3] for r in rs]
    out = {"n": len(rs), "mse": sum((p - y) ** 2 for p, y in zip(preds, ys)) / len(rs), "pearson": pearson(preds, ys)}
    # maternal - paternal difference per CpG x cross
    pair = {}
    for r, p in zip(rs, preds):
        pair.setdefault((r[0], r[1]), {})[r[2]] = (p, r[3], r[5], r[6])
    pd_, od_, loci = [], [], {}
    for k, v in pair.items():
        if 0 in v and 1 in v:
            pdiff = v[0][0] - v[1][0]; odiff = v[0][1] - v[1][1]
            pd_.append(pdiff); od_.append(odiff)
            if v[0][2]:
                loci.setdefault(v[0][3], [v[0][2], [], []])
                loci[v[0][3]][1].append(odiff); loci[v[0][3]][2].append(pdiff)
    out["diff_pearson"] = pearson(pd_, od_) if len(pd_) > 2 and torch.tensor(pd_).std() > 1e-6 else 0.0
    out["mean_abs_obs_diff"] = sum(abs(x) for x in od_) / max(len(od_), 1)
    out["mean_abs_pred_diff"] = sum(abs(x) for x in pd_) / max(len(pd_), 1)
    if loci:
        def sign_ok(par, m):
            return (m > 0) if par == "maternal" else (m < 0)
        obs_ok = sum(sign_ok(par, sum(o) / len(o)) for par, o, p in loci.values())
        pred_ok = sum(sign_ok(par, sum(p) / len(p)) for par, o, p in loci.values())
        out["dmr_loci"] = len(loci); out["dmr_obs_direction_matches_published"] = obs_ok; out["dmr_pred_direction_matches_published"] = pred_ok
    return out, preds

ymean = sum(r[3] for r in train) / len(train)
results = {"config": dict(win=WIN, train_pct=TRAIN_PCT, mincov=MINCOV, epochs=EPOCHS, batch=BATCH, n_train=len(train), n_test=len(test), n_dmr=len(dmr)),
           "baseline_predict_train_mean": {"test_mse": sum((ymean - r[3]) ** 2 for r in test) / len(test), "dmr_mse": sum((ymean - r[3]) ** 2 for r in dmr) / max(len(dmr), 1)}}
log("baseline", results["baseline_predict_train_mean"])

for cond, use_parent in (("A_sequence_only", False), ("B_sequence_plus_parent", True)):
    torch.manual_seed(0); random.seed(0)
    model = Net(use_parent).cuda()
    opt = torch.optim.AdamW([{"params": model.b.parameters(), "lr": 3e-5},
                             {"params": list(model.head.parameters()) + list(model.pe.parameters()), "lr": 1e-3}], weight_decay=0.01)
    steps = int(math.ceil(len(train) / BATCH) * EPOCHS)
    sched = torch.optim.lr_scheduler.OneCycleLR(opt, max_lr=[3e-5, 1e-3], total_steps=steps, pct_start=0.05)
    log(cond, "training", steps, "steps")
    t0 = time.time(); step = 0; run = 0.0
    while step < steps:
        for ids, am, par, y, _ in batches(train, BATCH, True):
            with torch.autocast("cuda", dtype=torch.bfloat16):
                loss = nn.functional.mse_loss(model(ids, am, par), y)
            opt.zero_grad(set_to_none=True); loss.backward(); opt.step(); sched.step()
            run = 0.98 * run + 0.02 * loss.item(); step += 1
            if step % 500 == 0:
                log(cond, "step", step, "loss(ema) %.4f" % run, "%.1f ex/s" % (step * BATCH / (time.time() - t0)))
            if step >= steps:
                break
    res = {"train_minutes": (time.time() - t0) / 60}
    res["test_chr19"], _ = evaluate(model, test)
    res["dmr"], _ = evaluate(model, dmr)
    results[cond] = res
    log(cond, json.dumps(res))
    del model, opt; torch.cuda.empty_cache()

json.dump(results, open(OUT + "/metrics.json", "w"), indent=2)
log("RESULTS", json.dumps(results, indent=2))
print("FINETUNE_OK")
