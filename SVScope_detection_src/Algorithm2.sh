#!/bin/bash

# ------------------------------
# Automate the Algorithm2's SV testing process
# ------------------------------

# default parameter
BAM=""
OUT_DIR=""
REFERENCE=""
JOBS=4  # Default number of threads

# Parsing command line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --bam_path) BAM="$2"; shift ;;
        --out_dir) OUT_DIR="$2"; shift ;;
        --referenceFasta) REFERENCE="$2"; shift ;;
        --jobs) JOBS="$2"; shift ;;
        *) echo "unknown parameter: $1" >&2; exit 1 ;;
    esac
    shift
done

# Checking parameters
if [[ -z "$BAM" || -z "$OUT_DIR" || -z "$REFERENCE" ]]; then
    echo "❌ The parameters are incomplete! Please provide：--bam_path --out_dir --referenceFasta [--jobs N]"
    exit 1
fi

# Checking configManta.py path
CONFIG_MANTA=$(which configManta.py)
if [[ -z "$CONFIG_MANTA" ]]; then
    echo "❌ Not found configManta.py，Make sure that you have activated the SVScope environment or set the PATH。"
    exit 1
fi
echo "🔧 Make sure that you have activated the SVScope environment or set the PATH：$CONFIG_MANTA"

# Creating an output directory
mkdir -p "$OUT_DIR"

echo "📁 Configuring SVScope Workflow..."
"$CONFIG_MANTA" \
    --bam "$BAM" \
    --exome \
    --runDir "$OUT_DIR" \
    --referenceFasta "$REFERENCE"

if [[ ! -f "$OUT_DIR/runWorkflow.py" ]]; then
    echo "❌ Configuration failed, not found runWorkflow.py！"
    exit 1
fi

echo "🚀 Run SVScope workflows using $JOBS threads..."
python "$OUT_DIR/runWorkflow.py" --jobs "$JOBS"

VCF_SRC="$OUT_DIR/results/variants/diploidSV.vcf.gz"
VCF_GUNZIPPED="$OUT_DIR/results/variants/diploidSV.vcf"
VCF_RENAMED="./Algorithm2.vcf"

# Checking and copying zip files
if [[ -f "$VCF_SRC" ]]; then
    cp "$VCF_SRC" "./$(basename "$VCF_SRC")"
    echo "✅ A VCF file has been generated and copied to the current directory：$(basename "$VCF_SRC")"

    # Unzip to the output directory
    gunzip -c "$VCF_SRC" > "$VCF_GUNZIPPED"

    if [[ -f "$VCF_GUNZIPPED" ]]; then
        # Then copy it as Algorithm2.vcf
        cp "$VCF_GUNZIPPED" "$VCF_RENAMED"
        echo "📄 Decompressed and copied to: $VCF_RENAMED"
    else
        echo "❌ Decompression failed: $VCF_GUNZIPPED not found"
        exit 1
    fi
else
    echo "❌ The resultant file diploidSV.vcf.gz was not found and SVScope may not have completed the run correctly."
    exit 1
fi