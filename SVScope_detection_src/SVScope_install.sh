#!/bin/bash

set -e  # Exit the script if there's an error
ENV_NAME="SVScope"
INSTALL_DIR="$HOME/SVScope_tools"
SURVIVOR_DIR="$INSTALL_DIR/SURVIVOR"

echo "🧬 [SVScope] Started installing a set of multi-source fusion algorithms for structural variation detection..."

# Step 1: Initializing conda shell support
if ! command -v conda >/dev/null 2>&1; then
    echo "❌ conda is not detected, please install Miniconda or Anaconda first!"
    exit 1
fi

source "$(conda info --base)/etc/profile.d/conda.sh"

# Step 2: Create the environment (if it does not already exist)
if conda info --envs | grep -q "$ENV_NAME"; then
    echo "🔄 Conda environment '$ENV_NAME' already exists, skip creation."
else
    echo "📦 The conda environment is being created. '$ENV_NAME'..."
    conda create -n "$ENV_NAME" -y
fi

# Step 3: activation environment
echo "🚀 The environment is being activated '$ENV_NAME'..."
conda activate "$ENV_NAME"

# Step 4: Installation of Algorithm 1 and Algorithm 2（Using bioconda）
echo "📥Algorithm 1 and Algorithm 2 are being installed......."
conda install -y -c bioconda manta

conda install -y -c bioconda lumpy-sv
conda install -y -c bioconda samblaster

# Step 5: Install Delly (manual download and link)
echo "📥 Algorithm 3 is being installed......."
mkdir -p "$INSTALL_DIR"
cd "$INSTALL_DIR"
wget -q https://github.com/dellytools/delly/releases/download/v1.0.3/delly_v1.0.3_linux_x86_64bit
chmod +x delly_v1.0.3_linux_x86_64bit
ln -sf "$INSTALL_DIR/delly_v1.0.3_linux_x86_64bit" "$INSTALL_DIR/delly"

# Step 6: Add Delly to PATH (current script only)
export PATH="$INSTALL_DIR:$PATH"


# Step 6.5: Permanently add delly to PATH
BASHRC="$HOME/.bashrc"
DELLY_PATH_LINE="export PATH=\"$INSTALL_DIR:\$PATH\""
if ! grep -Fxq "$DELLY_PATH_LINE" "$BASHRC"; then
    echo "" >> "$BASHRC"
    echo "# Add delly to PATH" >> "$BASHRC"
    echo "$DELLY_PATH_LINE" >> "$BASHRC"
    echo "✅ delly path has been added to your ~/.bashrc"
else
    echo "ℹ️ delly path is already in ~/.bashrc, skipping."
fi


# Step 7: Install SURVIVOR from source
echo "📥 Installing SURVIVOR..."
cd "$INSTALL_DIR"

if [ ! -d "$SURVIVOR_DIR" ]; then
    git clone https://github.com/fritzsedlazeck/SURVIVOR.git
fi

cd "$SURVIVOR_DIR/Debug"
make

# Step 8: Add SURVIVOR to PATH in current session
export PATH="$SURVIVOR_DIR/Debug:$PATH"
echo "✅ SURVIVOR temporarily added to PATH for current session."

# Step 9: Permanently add SURVIVOR to PATH by writing to .bashrc
# BASHRC="$HOME/.bashrc"
SURVIVOR_PATH_LINE="export PATH=\"$SURVIVOR_DIR/Debug:\$PATH\""

# Avoid duplicate entries
if ! grep -Fxq "$SURVIVOR_PATH_LINE" "$BASHRC"; then
    echo "" >> "$BASHRC"
    echo "# Add SURVIVOR to PATH" >> "$BASHRC"
    echo "$SURVIVOR_PATH_LINE" >> "$BASHRC"
    echo "✅ SURVIVOR path has been added to your ~/.bashrc"
else
    echo "ℹ️ SURVIVOR path is already in ~/.bashrc, skipping."
fi

# Optional: apply .bashrc changes immediately (works only in interactive shells)
source "$BASHRC"


# Step 9: Check if the tool was installed successfully
echo "✅ Installation complete, check tool availability："

echo -e "  ➤ Algorithm 1: \c"
command -v configManta.py >/dev/null 2>&1 && echo "Installed ✅" || echo "Not Found ❌"

echo -e "  ➤ Algorithm 2: \c"
command -v lumpyexpress >/dev/null 2>&1 && echo "Installed ✅" || echo "Not Found ❌"

echo -e "  ➤ Algorithm 3: \c"
command -v delly >/dev/null 2>&1 && echo "Installed ✅" || echo "Not Found ❌"

echo -e "  ➤ SURVIVOR: \c"
command -v SURVIVOR >/dev/null 2>&1 && echo "Installed ✅" || echo "Not Found ❌"

echo "🎉 [SVScope] All algorithm installation processes are finished! Make sure you activate the environment afterwards: conda activate $ENV_NAME"

echo “👉 Please run ‘source ~/.bashrc’ to validate the path manually to avoid interrupting the current conda environment”
