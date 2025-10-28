echo "# QSAR-Predictor

A command-line QSAR predictor using **RDKit** and **Mordred** descriptors for small molecules.  
Predicts **α-glucosidase inhibitory activity** (pIC50 regression + Active/Inactive classification).

---

## 🚀 Usage

Clone and install dependencies:
\`\`\`bash
git clone https://github.com/bryangervais/QSAR-Predictor.git
cd QSAR-Predictor
pip install -r requirements.txt
\`\`\`

Run prediction on a SMILES string:
\`\`\`bash
python predictor.py "COC1=CC=C(C=C1)O"
\`\`\`

Example output:
\`\`\`
SMILES: COC1=CC=C(C=C1)O
Predicted pIC50: 4.984
Predicted class: Inactive
\`\`\`

---

## 🧠 Model Details
- Descriptors: Mordred 2D descriptors + engineered phenolic_OH_count  
- Regressor: Random Forest Regressor  
- Classifier: Random Forest Classifier  
- Input: SMILES string  
- Output: Predicted pIC50 (continuous) and binary activity

---

## 📦 Files
- \`predictor.py\` — CLI tool  
- \`model/\` — trained models (`.pkl` files)  
- \`requirements.txt\` — dependencies  

---

## 🧑‍💻 Author
Created by **Bryan Gervais** ([@bryangervais](https://github.com/bryangervais))  
Version 1.0 — October 2025
" > README.md
git add README.md
git commit -m "docs: add README with usage instructions"
git push
