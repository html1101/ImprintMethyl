set -x
module load anaconda3
source activate base
conda create -y -q --clone progen -n dnabert2
source activate dnabert2
python -c "import torch,sys;print(\"torch\",torch.__version__,\"cuda\",torch.version.cuda,sys.version)"
pip install -q "transformers==4.29.2" einops sentencepiece huggingface_hub accelerate
export HF_HOME=${IMPRINT_BASE:-/ocean/projects/cis260266p/rpandit1/imprintmethyl}/hf_cache
python - <<PY
from huggingface_hub import snapshot_download
p=snapshot_download("zhihan1996/DNABERT-2-117M")
print("model at",p)
PY
echo SETUP_DONE
