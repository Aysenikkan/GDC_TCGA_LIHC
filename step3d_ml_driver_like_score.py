import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, roc_curve, classification_report, precision_recall_curve, average_precision_score

# ============================================================
# STEP 3D: Weak-supervised ML -> "driver-like" score
# Girdi : outputs/gene_priority_score.csv (STEP 2)
# Çıktı : outputs/step3d_ml_gene_scores.csv
#         outputs/step3d_roc_curve.png
#         outputs/step3d_pr_curve.png
#         outputs/step3d_report.txt
# ============================================================

BASE_DIR = r"D:\ALSU\GDC_TCGA_LIHC"   # gerekirse değiştir
OUTPUT_DIR = os.path.join(BASE_DIR, "outputs")
INPUT_PATH = os.path.join(OUTPUT_DIR, "gene_priority_score.csv")

OUT_SCORES = os.path.join(OUTPUT_DIR, "step3d_ml_gene_scores.csv")
OUT_ROC    = os.path.join(OUTPUT_DIR, "step3d_roc_curve.png")
OUT_PR     = os.path.join(OUTPUT_DIR, "step3d_pr_curve.png")
OUT_REPORT = os.path.join(OUTPUT_DIR, "step3d_report.txt")

os.makedirs(OUTPUT_DIR, exist_ok=True)

print("Okunan dosya:", INPUT_PATH)
df = pd.read_csv(INPUT_PATH)

print("Shape:", df.shape)
print("Kolonlar:", list(df.columns))

# ------------------------------------------------------------
# 1) Özellik kolonları (varsa hotspot_ratio kullan)
# ------------------------------------------------------------
base_features = ["patient_frequency", "high_impact_ratio", "n_mutations", "n_patients"]
if "hotspot_ratio" in df.columns:
    base_features.append("hotspot_ratio")

# log transform (mutasyon sayısı gene-length bias etkisini yumuşatır)
df["log_n_mutations"] = np.log1p(pd.to_numeric(df["n_mutations"], errors="coerce").fillna(0))
features = ["patient_frequency", "high_impact_ratio", "log_n_mutations", "n_patients"]
if "hotspot_ratio" in df.columns:
    features.append("hotspot_ratio")

# numeric'e çevir
for c in features:
    df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0)

# ------------------------------------------------------------
# 2) Weak label oluştur: mini driver list = 1, diğerleri = 0
#    (İstersen listeyi büyütürüz)
# ------------------------------------------------------------
known_driver_genes = set([
    "TP53", "CTNNB1", "AXIN1", "ARID1A", "ALB",
    "ACVR2A", "RB1", "RPS6KA3", "BAP1", "APOB",
    "CDKN2A", "KEAP1", "NFE2L2", "PIK3CA", "TERT"
])

df["weak_label_driver"] = df["Hugo_Symbol"].astype(str).isin(known_driver_genes).astype(int)

pos = int(df["weak_label_driver"].sum())
neg = int((df["weak_label_driver"] == 0).sum())
print(f"Weak labels -> Positive(driver): {pos}, Negative: {neg}")

# Çok az pozitif varsa öğrenme zorlaşır ama yine de çalışır.
# ------------------------------------------------------------
# 3) Train/Test split (stratify ile)
# ------------------------------------------------------------
X = df[features].values
y = df["weak_label_driver"].values

# Eğer pos çok azsa stratify hata verebilir; güvenli kontrol:
stratify = y if (pos >= 5) else None

X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.25, random_state=42, stratify=stratify
)

# ------------------------------------------------------------
# 4) Modeller
#    A) Logistic Regression (yorumlanabilir, hızlı)
#    B) Random Forest (non-linear)
# ------------------------------------------------------------
logreg = Pipeline([
    ("scaler", StandardScaler()),
    ("clf", LogisticRegression(max_iter=2000, class_weight="balanced", random_state=42))
])

rf = RandomForestClassifier(
    n_estimators=400,
    random_state=42,
    class_weight="balanced_subsample",
    max_depth=None
)

models = {
    "LogReg": logreg,
    "RandomForest": rf
}

results = {}

for name, model in models.items():
    model.fit(X_train, y_train)

    # probability
    if hasattr(model, "predict_proba"):
        prob_test = model.predict_proba(X_test)[:, 1]
    else:
        # fallback
        prob_test = model.decision_function(X_test)

    auc = roc_auc_score(y_test, prob_test) if len(np.unique(y_test)) > 1 else np.nan
    ap  = average_precision_score(y_test, prob_test) if len(np.unique(y_test)) > 1 else np.nan

    results[name] = {"model": model, "auc": auc, "ap": ap}

# En iyi modeli seç (AUC yüksek olan)
best_name = sorted(results.keys(), key=lambda k: (results[k]["auc"] if not np.isnan(results[k]["auc"]) else -1), reverse=True)[0]
best_model = results[best_name]["model"]

print("\nEn iyi model:", best_name, "AUC:", results[best_name]["auc"], "AP:", results[best_name]["ap"])

# ------------------------------------------------------------
# 5) Tüm genlere ML probability üret
# ------------------------------------------------------------
if hasattr(best_model, "predict_proba"):
    df["ml_driver_probability"] = best_model.predict_proba(X)[:, 1]
else:
    # normalize decision score to 0-1
    s = best_model.decision_function(X)
    df["ml_driver_probability"] = (s - s.min()) / (s.max() - s.min() + 1e-9)

# İstersen hibrit skor (eski score + ML)
# 0.6 ML + 0.4 gene_priority_score
if "gene_priority_score" in df.columns:
    df["hybrid_score"] = 0.6 * df["ml_driver_probability"] + 0.4 * pd.to_numeric(df["gene_priority_score"], errors="coerce").fillna(0)
else:
    df["hybrid_score"] = df["ml_driver_probability"]

# Sırala
df_sorted = df.sort_values("ml_driver_probability", ascending=False).reset_index(drop=True)

df_sorted.to_csv(OUT_SCORES, index=False)
print("\n✅ ML skor dosyası kaydedildi:")
print("->", OUT_SCORES)

# ------------------------------------------------------------
# 6) ROC + PR eğrileri (test set)
# ------------------------------------------------------------
# Test probability again
if hasattr(best_model, "predict_proba"):
    prob_test = best_model.predict_proba(X_test)[:, 1]
else:
    s = best_model.decision_function(X_test)
    prob_test = (s - s.min()) / (s.max() - s.min() + 1e-9)

# ROC
if len(np.unique(y_test)) > 1:
    fpr, tpr, _ = roc_curve(y_test, prob_test)
    plt.figure(figsize=(7, 5))
    plt.plot(fpr, tpr)
    plt.plot([0, 1], [0, 1], linestyle="--")
    plt.title(f"ROC Curve ({best_name})  AUC={results[best_name]['auc']:.3f}")
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.tight_layout()
    plt.savefig(OUT_ROC, dpi=200)
    plt.show()

    # PR
    prec, rec, _ = precision_recall_curve(y_test, prob_test)
    ap = average_precision_score(y_test, prob_test)
    plt.figure(figsize=(7, 5))
    plt.plot(rec, prec)
    plt.title(f"Precision-Recall ({best_name})  AP={ap:.3f}")
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.tight_layout()
    plt.savefig(OUT_PR, dpi=200)
    plt.show()

    print("\n✅ ROC ve PR grafikleri kaydedildi:")
    print("->", OUT_ROC)
    print("->", OUT_PR)
else:
    print("\n⚠️ Test set tek sınıf içeriyor, ROC/PR çizilemedi.")

# ------------------------------------------------------------
# 7) Mini rapor üret
# ------------------------------------------------------------
top20 = df_sorted[["Hugo_Symbol", "ml_driver_probability", "hybrid_score", "patient_frequency", "high_impact_ratio"] + (["hotspot_ratio"] if "hotspot_ratio" in df.columns else [])].head(20)

# known drivers top kaçta?
driver_ranks = []
sym_list = df_sorted["Hugo_Symbol"].astype(str).tolist()
for g in sorted(list(known_driver_genes)):
    if g in sym_list:
        driver_ranks.append((g, sym_list.index(g) + 1, float(df_sorted.loc[sym_list.index(g), "ml_driver_probability"])))

with open(OUT_REPORT, "w", encoding="utf-8") as f:
    f.write("STEP 3D REPORT - Weak-supervised ML driver-like scoring\n")
    f.write("======================================================\n\n")
    f.write(f"Input file: {INPUT_PATH}\n")
    f.write(f"Best model: {best_name}\n")
    f.write(f"AUC: {results[best_name]['auc']}\n")
    f.write(f"AP : {results[best_name]['ap']}\n\n")
    f.write(f"Weak labels -> Positive(driver): {pos}, Negative: {neg}\n")
    f.write(f"Features used: {features}\n\n")

    f.write("Top 20 genes by ML probability:\n")
    f.write("--------------------------------\n")
    f.write(top20.to_string(index=False))
    f.write("\n\n")

    f.write("Known driver mini-list ranks (gene, rank, prob):\n")
    f.write("-----------------------------------------------\n")
    for g, r, p in driver_ranks:
        f.write(f"- {g:10s} rank={r:4d}  prob={p:.4f}\n")

    f.write("\nNotes:\n")
    f.write("- This is weak-supervised learning (not true ground truth).\n")
    f.write("- Next improvement: use external curated driver lists (IntOGen, COSMIC Cancer Gene Census, OncoKB) as labels.\n")
    f.write("- Or do pan-cancer training then test on LIHC.\n")

print("\n✅ STEP 3D raporu kaydedildi:")
print("->", OUT_REPORT)

print("\nBİTTİ ✅ STEP 3D tamam.")
