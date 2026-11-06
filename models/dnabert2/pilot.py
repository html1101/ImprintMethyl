import os, sys, time, random, torch
os.environ.setdefault("HF_HOME", os.environ.get("IMPRINT_BASE", "/ocean/projects/cis260266p/rpandit1/imprintmethyl") + "/hf_cache")
NO_TRITON = os.environ.get("NO_TRITON")=="1"
if NO_TRITON:
    sys.modules["triton"] = None   # force DNABERT-2 to use its plain-PyTorch attention
from transformers import AutoTokenizer, AutoModel
name="zhihan1996/DNABERT-2-117M"
tok=AutoTokenizer.from_pretrained(name, trust_remote_code=True)
model=AutoModel.from_pretrained(name, trust_remote_code=True).cuda()
print("params(M):", sum(p.numel() for p in model.parameters())/1e6, "no_triton:", NO_TRITON, flush=True)
random.seed(0)
L=int(os.environ.get("SEQLEN","1000")); B=int(os.environ.get("BATCH","32"))
seqs=["".join(random.choice("ACGT") for _ in range(L)) for _ in range(B)]
enc=tok(seqs, return_tensors="pt", padding=True)
ids=enc["input_ids"].cuda(); am=enc["attention_mask"].cuda()
print("tokens per seq:", ids.shape[1], "=> bases/token:", round(L/ids.shape[1],2), flush=True)
head=torch.nn.Linear(768,1).cuda()
opt=torch.optim.AdamW(list(model.parameters())+list(head.parameters()), lr=2e-5)
y=torch.rand(B,1,device="cuda")
def step():
    with torch.autocast("cuda", dtype=torch.bfloat16):
        h=model(ids, attention_mask=am)[0]            # [B,T,768]
        pooled=(h*am.unsqueeze(-1)).sum(1)/am.sum(1,keepdim=True)
        loss=torch.nn.functional.mse_loss(torch.sigmoid(head(pooled.float())), y)
    opt.zero_grad(); loss.backward(); opt.step(); return loss.item()
print("first loss:", step(), flush=True)
for _ in range(3): step()
torch.cuda.synchronize(); n=20; t=time.time()
for _ in range(n): step()
torch.cuda.synchronize(); dt=(time.time()-t)/n
print(f"train step: {dt*1000:.1f} ms for batch {B} x {L} bp => {B/dt:.1f} examples/s; peak mem {torch.cuda.max_memory_allocated()/1e9:.1f} GB", flush=True)
print("PILOT_OK")
