import sys, json, warnings, argparse
import joblib, pandas as pd, numpy as np
from pathlib import Path
from rdkit import Chem
from mordred import Calculator, descriptors

MODEL_DIR = Path("model")
REG_PKL = MODEL_DIR / "AGI_all_rf_regressor.pkl"
CLS_PKL = MODEL_DIR / "AGI_all_rf_classifier.pkl"
FEATURES_JSON = MODEL_DIR / "features.json"   # optional

def compute_engineered_features(mol: Chem.Mol) -> dict:
    phenol_smarts = Chem.MolFromSmarts("[OX2H][c]")
    n = len(mol.GetSubstructMatches(phenol_smarts)) if (mol and phenol_smarts) else 0
    return {"phenolic_OH_count": float(n)}

def smiles_to_descriptors(smiles: str) -> pd.DataFrame:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    calc = Calculator(descriptors, ignore_3D=True)
    df = calc.pandas([mol]).fillna(0)
    if isinstance(df.columns, pd.MultiIndex):
        df.columns = ["_".join([str(c) for c in col if c != ""]) for col in df.columns]
    for k, v in compute_engineered_features(mol).items():
        df[k] = v
    return df

def get_expected_features(model) -> list[str] | None:
    names = getattr(model, "feature_names_in_", None)
    if names is not None:
        return list(names)
    if FEATURES_JSON.exists():
        with open(FEATURES_JSON, "r") as f:
            obj = json.load(f)
        return obj["features"] if isinstance(obj, dict) and "features" in obj else obj
    return None

def align_columns(df: pd.DataFrame, expected: list[str]) -> pd.DataFrame:
    missing = [c for c in expected if c not in df.columns]
    if missing:
        warnings.warn("Missing engineered/features: " + ", ".join(missing[:10]) + ("..." if len(missing) > 10 else ""))
        for c in missing:
            df[c] = 0.0
    return df.reindex(columns=expected, fill_value=0.0)

def load_models():
    return joblib.load(REG_PKL), joblib.load(CLS_PKL)

def predict_df(df: pd.DataFrame) -> pd.DataFrame:
    reg, clf = load_models()
    expected = get_expected_features(reg) or get_expected_features(clf)
    if expected is None:
        raise RuntimeError("Cannot determine training feature list; provide model.feature_names_in_ or model/features.json")
    X = align_columns(df, expected)
    y_reg = reg.predict(X)
    y_cls = clf.predict(X)
    out = pd.DataFrame({
        "Pred_pIC50": y_reg.astype(float),
        "Pred_class": np.where(y_cls.astype(int) == 1, "Active", "Inactive"),
    }, index=df.index)
    return out

def main():
    ap = argparse.ArgumentParser(description="QSAR predictor (single SMILES or CSV).")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--smiles", help="Single SMILES string")
    g.add_argument("--csv", help="CSV path with a column named 'smiles'")
    ap.add_argument("--out", help="Output CSV path (for batch mode)")
    args = ap.parse_args()

    if args.smiles:
        df = smiles_to_descriptors(args.smiles)
        res = predict_df(df)
        print(f"\nSMILES: {args.smiles}")
        print(f"Predicted pIC50: {float(res.iloc[0,0]):.3f}")
        print(f"Predicted class: {res.iloc[0,1]}")
    else:
        inpath = Path(args.csv)
        raw = pd.read_csv(inpath)
        if "smiles" not in raw.columns:
            raise ValueError("CSV must have a column named 'smiles'")
        feats = []
        idx = []
        for i, s in enumerate(raw["smiles"].astype(str)):
            try:
                feats.append(smiles_to_descriptors(s))
                idx.append(i)
            except Exception as e:
                feats.append(pd.DataFrame())  # placeholder
                idx.append(i)
                print(f"[WARN] row {i}: {e}")
        DF = pd.concat(feats, axis=0, ignore_index=True)
        DF.index = idx
        res = predict_df(DF)
        outpath = Path(args.out) if args.out else inpath.with_suffix(".pred.csv")
        res.to_csv(outpath, index=False)
        print(f"Saved predictions to: {outpath}")

if __name__ == "__main__":
    main()
