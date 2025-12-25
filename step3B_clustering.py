import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

# ============================================================
# STEP 3B: Unsupervised ML (KMeans) + Elbow + Silhouette
# Girdi : outputs/gene_priority_score.csv  (veya gene_feature_table.csv)
# Çıktı : outputs/step3b_kmeans_genes.csv + grafikler + kısa rapor
# ============================================================

BASE_DIR = r"D:\ALSU\GDC_TCGA_LIHC"   # <-- kendi yolun
INPUT_PATH = os.path.join(BASE_DIR, "outputs", "gene_priority_score.csv")
OUTPUT_DIR = os.path.join(BASE_DIR, "outputs")
os.makedirs(OUTPUT_DIR, exist_ok=True)

OUT_CSV = os.path.join(OUTPUT_DIR, "step3b_kmeans_genes.csv")
PLOT_ELBOW = os.path.join(OUTPUT_DIR, "step3b_elbow_inertia.png")
PLOT_SIL = os.path.join(OUTPUT_DIR, "step3b_silhouette_scores.png")
REPORT_TXT = os.path.join(OUTPUT_DIR, "step3b_report.txt")

print("Okunan dosya:", INPUT_PATH)
df = pd.read_csv(INPUT_PATH)
print("Shape:", df.shape)
print("Kolonlar:", list(df.columns))

# ------------------------------------------------------------
# 1) Cluster feature seti seçelim
# (skorlar + oranlar; gen ismi hariç)
# ------------------------------------------------------------
feature_cols = [
    "n_mutations",
    "n_patients",
    "high_impact_ratio",
    "patient_frequency",
    "gene_priority_score"
]

# Eğer hotspot_ratio varsa ekleyelim
if "hotspot_ratio" in df.columns:
    feature_cols.append("hotspot_ratio")

# güvenlik
missing = [c for c in feature_cols if c not in df.columns]
if missing:
    raise ValueError(f"Feature kolonları eksik: {missing}\nMevcut kolonlar: {list(df.columns)}")

X = df[feature_cols].copy()

# numeric'e zorlama
for c in feature_cols:
    X[c] = pd.to_numeric(X[c], errors="coerce")
X = X.fillna(0)

# ------------------------------------------------------------
# 2) Ölçekleme (çok önemli!)
# ------------------------------------------------------------
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# ------------------------------------------------------------
# 3) Elbow (Inertia) + Silhouette hesapla
# ------------------------------------------------------------
k_min, k_max = 2, 12
k_values = list(range(k_min, k_max + 1))

inertias = []
sil_scores = []

for k in k_values:
    km = KMeans(n_clusters=k, random_state=42, n_init=10)
    labels = km.fit_predict(X_scaled)
    inertias.append(km.inertia_)
    sil_scores.append(silhouette_score(X_scaled, labels))

# ------------------------------------------------------------
# 4) Basit "dirsek/knee" bulma (çizgiye en uzak nokta)
#   - inertia grafiğinde en iyi elbow tahmini
# ------------------------------------------------------------
def find_knee_point(x, y):
    """
    x: k değerleri
    y: inertia değerleri
    Çizginin (ilk nokta - son nokta) en uzağındaki noktayı knee kabul eder.
    """
    x = np.array(x)
    y = np.array(y)

    # iki uç nokta
    x1, y1 = x[0], y[0]
    x2, y2 = x[-1], y[-1]

    # her noktanın doğruya uzaklığı
    # |(y2-y1)x - (x2-x1)y + x2*y1 - y2*x1| / sqrt((y2-y1)^2 + (x2-x1)^2)
    num = np.abs((y2 - y1) * x - (x2 - x1) * y + (x2 * y1 - y2 * x1))
    den = np.sqrt((y2 - y1)**2 + (x2 - x1)**2)
    dist = num / den

    knee_index = int(np.argmax(dist))
    return int(x[knee_index])

k_elbow = find_knee_point(k_values, inertias)

# Elbow çevresinde silhouette ile ince ayar: (k_elbow-2 .. k_elbow+2) aralığında en iyi silhouette
candidate_window = [k for k in k_values if (k_elbow - 2) <= k <= (k_elbow + 2)]
best_k = max(candidate_window, key=lambda k: sil_scores[k_values.index(k)])

print("\nElbow (knee) tahmini:", k_elbow)
print("Elbow çevresinde en iyi silhouette veren k:", best_k)

# ------------------------------------------------------------
# 5) Grafikler (outputs'a kaydet)
# ------------------------------------------------------------
plt.figure(figsize=(8, 5))
plt.plot(k_values, inertias, marker="o")
plt.axvline(k_elbow, linestyle="--")
plt.title("Elbow Method (Inertia / WCSS)")
plt.xlabel("k")
plt.ylabel("Inertia (WCSS)")
plt.tight_layout()
plt.savefig(PLOT_ELBOW, dpi=200)
plt.show()

plt.figure(figsize=(8, 5))
plt.plot(k_values, sil_scores, marker="o")
plt.axvline(best_k, linestyle="--")
plt.title("Silhouette Scores by k")
plt.xlabel("k")
plt.ylabel("Silhouette Score")
plt.tight_layout()
plt.savefig(PLOT_SIL, dpi=200)
plt.show()

print("\n✅ Grafikler kaydedildi:")
print(" -", PLOT_ELBOW)
print(" -", PLOT_SIL)

# ------------------------------------------------------------
# 6) Final KMeans (best_k ile)
# ------------------------------------------------------------
final_km = KMeans(n_clusters=best_k, random_state=42, n_init=10)
df["cluster"] = final_km.fit_predict(X_scaled)

# cluster özet
cluster_counts = df["cluster"].value_counts().sort_index()

# cluster bazlı feature ortalamaları
cluster_means = df.groupby("cluster")[feature_cols].mean()

# kaydet
df.to_csv(OUT_CSV, index=False)

# rapor yaz
with open(REPORT_TXT, "w", encoding="utf-8") as f:
    f.write("STEP 3B - KMeans Clustering (Genes)\n")
    f.write("=================================\n\n")
    f.write(f"k range: {k_min}-{k_max}\n")
    f.write(f"Elbow(knee) suggested k: {k_elbow}\n")
    f.write(f"Final chosen k (elbow-window + best silhouette): {best_k}\n\n")

    f.write("Silhouette scores:\n")
    for k, s in zip(k_values, sil_scores):
        f.write(f"- k={k:2d}  silhouette={s:.4f}\n")

    f.write("\nInertias:\n")
    for k, inn in zip(k_values, inertias):
        f.write(f"- k={k:2d}  inertia={inn:.2f}\n")

    f.write("\nCluster counts:\n")
    for c, cnt in cluster_counts.items():
        f.write(f"- cluster {c}: {cnt} genes\n")

    f.write("\nCluster means (feature averages):\n")
    f.write(cluster_means.to_string())
    f.write("\n")

print("\n✅ STEP 3B tamamlandı.")
print("-> Cluster'lı çıktı:", OUT_CSV)
print("-> Rapor:", REPORT_TXT)
print("\nCluster counts:\n", cluster_counts)

# ============================================================
# 7) EK: Zengin özetler + ekstra grafikler (outputs'a kaydet)
# ============================================================

# Ek çıktı yolları
PLOT_CLUSTER_SIZES = os.path.join(OUTPUT_DIR, "step3b_cluster_sizes.png")
PLOT_PCA = os.path.join(OUTPUT_DIR, "step3b_pca_clusters.png")
TOP20_CSV = os.path.join(OUTPUT_DIR, "step3b_top20_with_clusters.csv")

# --- Top 20 gen + cluster bilgisi
top20 = df.sort_values("gene_priority_score", ascending=False).head(20).copy()
top20.to_csv(TOP20_CSV, index=False)

print("\n📌 Top 20 gen + cluster kaydedildi:")
print("->", TOP20_CSV)
print(top20[["Hugo_Symbol", "gene_priority_score", "cluster"]].head(20))

# --- Cluster boyutları grafiği
cluster_counts = df["cluster"].value_counts().sort_index()

plt.figure(figsize=(8, 5))
plt.bar(cluster_counts.index.astype(str), cluster_counts.values)
plt.title("Cluster Sizes (Number of Genes)")
plt.xlabel("Cluster")
plt.ylabel("Gene count")
plt.tight_layout()
plt.savefig(PLOT_CLUSTER_SIZES, dpi=200)
plt.show()

print("\n✅ Cluster size grafiği kaydedildi:")
print("->", PLOT_CLUSTER_SIZES)

# --- PCA ile 2D görselleştirme
from sklearn.decomposition import PCA

pca = PCA(n_components=2, random_state=42)
X_2d = pca.fit_transform(X_scaled)

plt.figure(figsize=(8, 6))
plt.scatter(X_2d[:, 0], X_2d[:, 1], s=10, alpha=0.6, c=df["cluster"])
plt.title("PCA (2D) - Genes colored by cluster")
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.tight_layout()
plt.savefig(PLOT_PCA, dpi=200)
plt.show()

print("\n✅ PCA grafiği kaydedildi:")
print("->", PLOT_PCA)

# --- Cluster bazlı özet istatistikleri (mean + median)
cluster_summary_mean = df.groupby("cluster")[feature_cols].mean()
cluster_summary_median = df.groupby("cluster")[feature_cols].median()

print("\nCluster means:\n", cluster_summary_mean)
print("\nCluster medians:\n", cluster_summary_median)

# --- Rapor dosyasına ekleme
with open(REPORT_TXT, "a", encoding="utf-8") as f:
    f.write("\n\nEXTRA OUTPUTS\n")
    f.write("=============\n")
    f.write(f"Top20 with clusters: {TOP20_CSV}\n")
    f.write(f"Cluster sizes plot:  {PLOT_CLUSTER_SIZES}\n")
    f.write(f"PCA plot:            {PLOT_PCA}\n\n")

    f.write("Cluster summary (MEAN):\n")
    f.write(cluster_summary_mean.to_string())
    f.write("\n\nCluster summary (MEDIAN):\n")
    f.write(cluster_summary_median.to_string())
    f.write("\n")
