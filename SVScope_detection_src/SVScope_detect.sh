#!/bin/bash

# ========== parameter analysis ==========
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam_path)
      BAM="$2"
      shift 2
      ;;
    --referenceFasta)
      REFERENCE="$2"
      shift 2
      ;;
    --out_dir)
      OUTDIR="$2"
      shift 2
      ;;
    --prefix)
      PREFIX="$2"
      shift 2
      ;;
    --jobs)
      JOBS="$2"
      shift 2
      ;;
    *)
      echo "❌ Unknown parameter: $1"
      exit 1
      ;;
  esac
done

# ========== Parameter checking ==========
if [[ -z "$BAM" || -z "$REFERENCE" || -z "$OUTDIR" || -z "$PREFIX" ]]; then
  echo "❗ Missing parameters. Required: --bam_path --referenceFasta --out_dir --prefix"
  exit 1
fi

mkdir -p "$OUTDIR"

echo "🚀 Launching SV detection with 3 algorithms in parallel..."

# ==========Parallel operation algorithm==========
bash Algorithm1.sh --bam_path "$BAM" --referenceFasta "$REFERENCE" --out_dir "$OUTDIR" &
PID1=$!

bash Algorithm2.sh --bam_path "$BAM" --referenceFasta "$REFERENCE" --out_dir "$OUTDIR" --jobs "${JOBS:-4}" &
PID2=$!

bash Algorithm3.sh --bam_path "$BAM" --out_dir "$OUTDIR" --prefix "$PREFIX" &
PID3=$!

wait $PID1 $PID2 $PID3

echo "✅ All algorithms completed. Proceeding to VCF sort and merge..."

# ========== Sorting VCF files==========
cd "$OUTDIR"

for i in 1 2 3; do
  VCF="Algorithm${i}.vcf"
  SORTED="Algorithm${i}_sorted.vcf"
  if [[ -f "$VCF" ]]; then
    echo "🔃 Sorting $VCF ..."
    bcftools sort "$VCF" -o "$SORTED"
    echo "✅ Sorted file generated: $SORTED"
  else
    echo "❌ Missing $VCF. Skipping..."
  fi
done

# ========== Build sample_file==========
echo "📄 Creating sample_file for SURVIVOR..."
SAMPLE_LIST="sample_file"
> "$SAMPLE_LIST"
for i in 1 2 3; do
  SORTED="Algorithm${i}_sorted.vcf"
  [[ -f "$SORTED" ]] && echo "$OUTDIR/$SORTED" >> "$SAMPLE_LIST"
done

# ========== Checking the SURVIVOR command==========
if ! command -v SURVIVOR &> /dev/null; then
  echo "❌ SURVIVOR is not in PATH. Please ensure the following line is in your ~/.bashrc and re-source it:"
  echo "export PATH=\"/full/path/to/SURVIVOR/Debug:\$PATH\""
  exit 1
fi

# ========== Run SURVIVOR merge==========
echo "🔗 Running SURVIVOR merge..."
SURVIVOR_OUTPUT="sample_merge.vcf"
SURVIVOR merge "$SAMPLE_LIST" 1000 1 1 1 0 30 "$SURVIVOR_OUTPUT"

if [[ -f "$SURVIVOR_OUTPUT" ]]; then
  echo "✅ SURVIVOR merged VCF generated: $SURVIVOR_OUTPUT"

  echo "🔃 Sorting final merged VCF..."
  bcftools sort "$SURVIVOR_OUTPUT" -o SVScope_unfiltered.vcf

  if [[ -f SVScope_unfiltered.vcf ]]; then
    echo "🎉 Final VCF is ready: $OUTDIR/SVScope_unfiltered.vcf"
  else
    echo "❌ Sorting failed on merged file."
  fi
else
  echo "❌ SURVIVOR merge failed!"
fi
