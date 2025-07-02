#!/bin/bash

# Set the default parameter to null
BAM=""
REFERENCE=""
OUTDIR=""

# parameterization
while [[ $# -gt 0 ]]; do
  case $1 in
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
    *)
      echo "❌ Unknown parameter: $1"
      exit 1
      ;;
  esac
done

# Check if the parameter is empty
if [[ -z "$BAM" || -z "$REFERENCE" || -z "$OUTDIR" ]]; then
  echo "❗ Required parameters are missing. Usage examples："
  echo "bash detection_delly.sh --bam_path sample.bam --referenceFasta ref.fa --out_dir ./output"
  exit 1
fi

# Creating an output directory
mkdir -p "$OUTDIR"

# Remove trailing slash from OUTDIR if exists
OUTDIR="${OUTDIR%/}"

# Get file prefix
PREFIX=$(basename "$BAM" .bam)
OUTBCF="$OUTDIR/${PREFIX}.delly.bcf"
OUTVCF="$OUTDIR/Algorithm1.vcf"

echo "🧬 BAM file：$BAM"
echo "🧬 Reference genome：$REFERENCE"
echo "📁 Output directory：$OUTDIR"
echo "📤 Output BCF：$OUTBCF"

# Check BAM for index files
if [[ ! -f "${BAM}.bai" ]]; then
  echo "🔍 BAM index not detected, being generated automatically..."
  samtools index "$BAM"
fi

# Run DELLY
echo "🚀 Running SVScope SV detection..."
delly call -o "$OUTBCF" -g "$REFERENCE" "$BAM"


# Check if the BCF was successfully generated
if [[ -f "$OUTBCF" ]]; then
  echo "✅ The BCF file was generated successfully：$OUTBCF"
  
  # Convert to VCF and rename to Algorithm1.vcf
  echo "🔄 Converting to VCF format..."
  bcftools view "$OUTBCF" > "$OUTVCF"
  
  if [[ -f "$OUTVCF" ]]; then
    echo "✅ VCF file generated successfully：$OUTVCF"
  else
    echo "❌ Conversion to VCF failed！"
  fi
else
  echo "❌ BCF file not successfully generated, SVScope detection failed！"
fi
