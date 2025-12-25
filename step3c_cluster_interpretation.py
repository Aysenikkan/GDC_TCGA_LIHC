import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# STEP 3C: Cluster interpretation + auto-labeling + mini report
# Girdi : outputs/step3b_kmeans_genes.csv
# Çıktı : outputs/step3c_cluster_summary.csv
#         outputs/step3c_cluster_labels.csv
#         outputs/step3c_cluster_interpretation_report.txt
#         outputs/step3c_score_by_cluster.png
# ============================================================

BASE_DIR = r"D:\ALSU\GDC_TCGA_LIHC"   # <-- gerekirse değiştir
OUTPUT_DIR = os.path.join(BASE_DIR, "outputs")
INPUT_PATH = os.path.join(OUTPUT_DIR, "step3b_kmeans_genes.csv")

SUMMARY_CSV = os.path.join(OUTPUT_DIR, "step3c_cluster_summary.csv")
LABELS_CSV  = os.path.join(OUTPUT_DIR, "step3c_cluster_labels.csv")
REPORT_TXT  = os.path.join(OUTPUT_DIR, "step3c_cluster_interpretation_report.txt")
PLOT_SCORE  = os.path.join(OUTPUT_DIR, "step3c_score_by_cluster.png")

os.makedirs(OUTPUT_DIR, exist_ok=True)

print("Okunan dosya:", INPUT_PATH)
df = pd.read_csv(INPUT_PATH)

print("Shape:", df.shape)
print("Kolonlar:", list(df.columns))

# ------------------------------------------------------------
# 1) Gerekli kolon kontrolü
# ------------------------------------------------------------
required = [
    "Hugo_Symbol",
    "cluster",
    "n_mutations",
    "n_patients",
    "patient_frequency",
    "high_impact_ratio",
    "gene_priority_score"
]
missing = [c for c in required if c not in df.columns]
if missing:
    raise ValueError(f"Eksik kolonlar var: {missing}\nMevcut kolonlar: {list(df.columns)}")

# Hotspot_ratio yoksa sorun değil; varsa kullanacağız
has_hotspot_ratio = "hotspot_ratio" in df.columns

# Sayısal kolonları düzelt
num_cols = ["n_mutations", "n_patients", "patient_frequency", "high_impact_ratio", "gene_priority_score"]
if has_hotspot_ratio:
    num_cols.append("hotspot_ratio")

for c in num_cols:
    df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0)

df["cluster"] = pd.to_numeric(df["cluster"], errors="coerce").fillna(-1).astype(int)

# ------------------------------------------------------------
# 2) Cluster özetleri (mean/median, size, top genes)
# ------------------------------------------------------------
feature_cols = ["n_mutations", "n_patients", "patient_frequency", "high_impact_ratio", "gene_priority_score"]
if has_hotspot_ratio:
    feature_cols.append("hotspot_ratio")

cluster_size = df.groupby("cluster").size().rename("n_genes")

cluster_mean = df.groupby("cluster")[feature_cols].mean().add_prefix("mean_")
cluster_median = df.groupby("cluster")[feature_cols].median().add_prefix("median_")

summary = pd.concat([cluster_size, cluster_mean, cluster_median], axis=1).reset_index()
summary.to_csv(SUMMARY_CSV, index=False)

print("\n✅ Cluster summary kaydedildi:")
print("->", SUMMARY_CSV)
print(summary.head())

# ------------------------------------------------------------
# 3) Otomatik cluster etiketi (heuristic, basit ve sunumluk)
# Mantık:
# - patient_frequency yüksekse: "Common across patients"
# - high_impact_ratio yüksekse: "High impact"
# - n_mutations çok yüksekse: "Many mutations (possible gene-length bias)"
# - (varsa) hotspot_ratio yüksekse: "Hotspot-enriched"
# ------------------------------------------------------------

# Normalize için yardımcı
def minmax(s):
    s = s.astype(float)
    mn, mx = s.min(), s.max()
    if mx - mn == 0:
        return pd.Series(np.zeros(len(s)), index=s.index)
    return (s - mn) / (mx - mn)

# summary üzerinden karar verelim
tmp = summary.copy()
tmp["score_common"] = minmax(tmp["mean_patient_frequency"])
tmp["score_impact"] = minmax(tmp["mean_high_impact_ratio"])
tmp["score_mut"]    = minmax(tmp["mean_n_mutations"])

if has_hotspot_ratio:
    tmp["score_hotspot"] = minmax(tmp["mean_hotspot_ratio"])
else:
    tmp["score_hotspot"] = 0.0

labels = []
reasons = []

for _, row in tmp.iterrows():
    c = int(row["cluster"])

    # En baskın özellik hangisi?
    scores = {
        "common": row["score_common"],
        "impact": row["score_impact"],
        "mut": row["score_mut"],
        "hotspot": row["score_hotspot"]
    }
    best = max(scores, key=scores.get)

    # Etiket üret
    if best == "common":
        label = "Common across patients"
        reason = "mean patient_frequency is relatively high"
    elif best == "impact":
        label = "High-impact enriched"
        reason = "mean high_impact_ratio is relatively high"
    elif best == "hotspot":
        label = "Hotspot-enriched"
        reason = "mean hotspot_ratio is relatively high"
    else:
        label = "Many mutations (check gene-length bias)"
        reason = "mean n_mutations is relatively high"

    labels.append(label)
    reasons.append(reason)

label_df = pd.DataFrame({
    "cluster": tmp["cluster"].astype(int),
    "auto_label": labels,
    "reason": reasons,
    "n_genes": tmp["n_genes"].astype(int),
    "mean_gene_priority_score": tmp["mean_gene_priority_score"].round(4),
    "mean_patient_frequency": tmp["mean_patient_frequency"].round(4),
    "mean_high_impact_ratio": tmp["mean_high_impact_ratio"].round(4),
    "mean_n_mutations": tmp["mean_n_mutations"].round(2),
})

if has_hotspot_ratio:
    label_df["mean_hotspot_ratio"] = tmp["mean_hotspot_ratio"].round(4)

label_df = label_df.sort_values("mean_gene_priority_score", ascending=False).reset_index(drop=True)
label_df.to_csv(LABELS_CSV, index=False)

print("\n✅ Cluster labels kaydedildi:")
print("->", LABELS_CSV)
print(label_df)

# ------------------------------------------------------------
# 4) Gen tablosuna label ekle + cluster içi top genler
# ------------------------------------------------------------
df2 = df.merge(label_df[["cluster", "auto_label"]], on="cluster", how="left")

# Her cluster için top 10 gene (score'a göre)
top_per_cluster = (
    df2.sort_values(["cluster", "gene_priority_score"], ascending=[True, False])
       .groupby("cluster")
       .head(10)
       [["cluster", "auto_label", "Hugo_Symbol", "gene_priority_score", "patient_frequency", "high_impact_ratio", "n_mutations"]]
)

# ------------------------------------------------------------
# 5) Grafik: Cluster bazlı score dağılımı (mean score)
# ------------------------------------------------------------
plot_df = df2.groupby("cluster")["gene_priority_score"].mean().sort_values(ascending=False)

plt.figure(figsize=(9, 5))
plt.bar(plot_df.index.astype(str), plot_df.values)
plt.title("Mean Gene Priority Score by Cluster")
plt.xlabel("Cluster")
plt.ylabel("Mean gene_priority_score")
plt.tight_layout()
plt.savefig(PLOT_SCORE, dpi=200)
plt.show()

print("\n✅ Grafik kaydedildi:")
print("->", PLOT_SCORE)

# ------------------------------------------------------------
# 6) Sunumluk rapor üret (TXT)
# ------------------------------------------------------------
with open(REPORT_TXT, "w", encoding="utf-8") as f:
    f.write("STEP 3C REPORT - LIHC Cluster Interpretation\n")
    f.write("==========================================\n\n")
    f.write(f"Input: {INPUT_PATH}\n")
    f.write(f"Total genes: {df2.shape[0]}\n")
    f.write(f"Total clusters: {df2['cluster'].nunique()}\n\n")

    f.write("Cluster labels (auto):\n")
    f.write("----------------------\n")
    f.write(label_df.to_string(index=False))
    f.write("\n\n")

    f.write("Top 10 genes per cluster (by gene_priority_score):\n")
    f.write("--------------------------------------------------\n")
    for c in sorted(df2["cluster"].unique()):
        sub = top_per_cluster[top_per_cluster["cluster"] == c]
        f.write(f"\nCluster {c}  |  {sub['auto_label'].iloc[0] if len(sub)>0 else 'NA'}\n")
        f.write(sub.to_string(index=False))
        f.write("\n")

    f.write("\nInterpretation guide:\n")
    f.write("- 'Common across patients' -> many patients carry mutations in these genes.\n")
    f.write("- 'High-impact enriched' -> higher ratio of high-impact variants.\n")
    if has_hotspot_ratio:
        f.write("- 'Hotspot-enriched' -> hotspot_ratio is relatively high.\n")
    f.write("- 'Many mutations' -> very mutation-heavy genes; check gene-length bias (e.g., TTN-like).\n")

print("\n✅ STEP 3C rapor dosyası oluşturuldu:")
print("->", REPORT_TXT)

print("\nBİTTİ ✅ 3C tamam.")
print("Sonraki adım: 3D (istersen) = ML ile gen skorunu iyileştirme / driver tahmin modeli")
