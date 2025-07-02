#!/bin/bash
set -euo pipefail

usage() {
    echo "Usage: $0 --bam_path BAM file path [--out_dir output directory] [--prefix (linguistics)]"
    echo "typical example: $0 --bam_path /path/to/sample.bam --out_dir /path/to/output --prefix sample"
    exit 1
}

# default value
OUTDIR="."
PREFIX="sample"

# parameter analysis
while [[ $# -gt 0 ]]; do
    case $1 in
        --bam_path)
            BAM="$2"
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
        *)
            echo "unknown parameter：$1"
            usage
            ;;
    esac
done

# Checking mandatory parameters
if [[ -z "${BAM:-}" ]]; then
    echo "Error: must be provided --bam_path parameters"
    usage
fi

# Make sure the output directory exists
mkdir -p "$OUTDIR"

echo "▶ Input BAM file：$BAM"
echo "▶ output directory：$OUTDIR"
echo "▶ File Prefix：$PREFIX"

# Tool checking
echo "🔍 Checking for necessary tools..."
TOOLS=("samtools" "lumpyexpress" "extractSplitReads_BwaMem")

for tool in "${TOOLS[@]}"; do
    TOOL_PATH=$(which "$tool" 2>/dev/null || true)
    if [[ -z "$TOOL_PATH" ]]; then
        echo "❌ Tool not found：$tool， Please make sure that you have installed and added the PATH"
        exit 1
    else
        echo "✅ find $tool ：$TOOL_PATH"
    fi
done

echo "🧩  Extraction  discordant paired-end reads..."
samtools view -b -F 1294 "$BAM" > "$OUTDIR/${PREFIX}.discordants.unsorted.bam" &

echo "🔪  Extraction split reads..."
samtools view -h "$BAM" \
    | extractSplitReads_BwaMem -i stdin \
    | samtools view -Sb - \
    > "$OUTDIR/${PREFIX}.splitters.unsorted.bam" &

# Wait for both background tasks to finish
wait

echo "🎉 Two files are generated to start sorting..."

samtools sort "$OUTDIR/${PREFIX}.discordants.unsorted.bam" -o "$OUTDIR/${PREFIX}.discordants.bam"
samtools sort "$OUTDIR/${PREFIX}.splitters.unsorted.bam" -o "$OUTDIR/${PREFIX}.splitters.bam"

# run lumpyexpress
echo "⚡ Run SVScope detecting structural variants..."
lumpyexpress \
    -B "$BAM" \
    -S "$OUTDIR/${PREFIX}.splitters.bam" \
    -D "$OUTDIR/${PREFIX}.discordants.bam" \
    -o "$OUTDIR/${PREFIX}.vcf"


VCF_PATH="$OUTDIR/${PREFIX}.vcf"
ALGO_VCF="$OUTDIR/Algorithm3.vcf"

if [[ -f "$VCF_PATH" ]]; then
  echo "✅ SV detection is complete. Result file: $VCF_PATH"
  
  # Copy and rename
  cp "$VCF_PATH" "$ALGO_VCF"
  echo "📄 Copied result as: $ALGO_VCF"
else
  echo "❌ SV detection failed or VCF file not found: $VCF_PATH"
fi
