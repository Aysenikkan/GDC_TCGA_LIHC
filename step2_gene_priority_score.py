import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# STEP 2: Gene priority score (mutasyon özelliklerinden skor)
# Girdi: outputs/gene_feature_table.csv
# Çıktı: outputs/gene_priority_score.csv + grafikler
# ============================================================

BASE_DIR = r"D:\ALSU\GDC_TCGA_LIHC"   # <-- KENDİ YOLUN FARKLIYSA DEĞİŞTİR

INPUT_PATH = os.path.join(BASE_DIR, "outputs", "gene_feature_table.csv")
OUTPUT_DIR = os.path.join(BASE_DIR, "outputs")
OUTPUT_PATH = os.path.join(OUTPUT_DIR, "gene_priority_score.csv")

os.makedirs(OUTPUT_DIR, exist_ok=True)

print("Okunan dosya:", INPUT_PATH)

# ------------------------------------------------------------
# 1) Gene feature tablosunu oku
# ------------------------------------------------------------
df = pd.read_csv(INPUT_PATH)

print("\nGene feature table shape:", df.shape)
print("Kolonlar:", list(df.columns))

# ------------------------------------------------------------
# 2) Zorunlu kolon kontrolü (hotspot_ratio dosyada yok, biz üreteceğiz)
# ------------------------------------------------------------
required_cols = ["Hugo_Symbol", "n_mutations", "n_patients", "hotspot_count", "high_impact_ratio", "patient_frequency"]
missing = [c for c in required_cols if c not in df.columns]

if missing:
    raise ValueError(f"Eksik kolonlar var: {missing}\n"
                     f"Mevcut kolonlar: {list(df.columns)}")

# ------------------------------------------------------------
# 3) Sayısal kolonları güvenli şekilde numeric'e çevir
# ------------------------------------------------------------
numeric_cols = ["n_mutations", "n_patients", "hotspot_count", "high_impact_ratio", "patient_frequency"]
for c in numeric_cols:
    df[c] = pd.to_numeric(df[c], errors="coerce")

df = df.dropna(subset=["Hugo_Symbol"])  # gen adı boşsa at
df = df.fillna(0)  # numeric NaN -> 0

# ------------------------------------------------------------
# 4) hotspot_ratio üret (bölme hatalarına karşı güvenli)
# ------------------------------------------------------------
# n_mutations 0 ise oranı 0 yap
df["hotspot_ratio"] = np.where(df["n_mutations"] > 0,
                               df["hotspot_count"] / df["n_mutations"],
                               0)

# ------------------------------------------------------------
# 5) Normalize (0-1 arası) fonksiyonu
# ------------------------------------------------------------
def minmax(series: pd.Series) -> pd.Series:
    s = series.astype(float)
    mn, mx = s.min(), s.max()
    if mx - mn == 0:
        return pd.Series(np.zeros(len(s)), index=s.index)
    return (s - mn) / (mx - mn)

# ------------------------------------------------------------
# 6) Skor metrikleri
#   0.50 * patient_frequency_norm
# + 0.30 * high_impact_ratio_norm
# + 0.20 * hotspot_ratio_norm
# ------------------------------------------------------------
df["patient_frequency_norm"] = minmax(df["patient_frequency"])
df["high_impact_ratio_norm"] = minmax(df["high_impact_ratio"])
df["hotspot_ratio_norm"] = minmax(df["hotspot_ratio"])

w_patient = 0.50
w_impact  = 0.30
w_hotspot = 0.20

df["gene_priority_score"] = (
    w_patient * df["patient_frequency_norm"] +
    w_impact  * df["high_impact_ratio_norm"] +
    w_hotspot * df["hotspot_ratio_norm"]
)

df["log_n_mutations"] = np.log1p(df["n_mutations"])

# ------------------------------------------------------------
# 7) Skora göre sırala ve kaydet
# ------------------------------------------------------------
df_sorted = df.sort_values("gene_priority_score", ascending=False).reset_index(drop=True)
df_sorted.to_csv(OUTPUT_PATH, index=False)

print("\n✅ Gene priority score oluşturuldu ve kaydedildi:")
print("->", OUTPUT_PATH)

print("\nToplam gen sayısı:", df_sorted.shape[0])

print("\n📌 Top 20 gen (skora göre):")
print(df_sorted[["Hugo_Symbol", "gene_priority_score", "patient_frequency", "high_impact_ratio", "hotspot_ratio",
                 "n_mutations", "n_patients"]].head(20))

# ------------------------------------------------------------
# 8) Grafikler (EKRANA GÖSTER + outputs klasörüne KAYDET)
# ------------------------------------------------------------

# Grafik 1: Top 20 gen - barplot
top20 = df_sorted.head(20).copy()

plt.figure(figsize=(12, 6))
plt.bar(top20["Hugo_Symbol"], top20["gene_priority_score"])
plt.xticks(rotation=75, ha="right")
plt.title("Top 20 Gene Priority Score (LIHC)")
plt.xlabel("Gene")
plt.ylabel("Gene Priority Score")
plt.tight_layout()

# Kaydet
plot1_path = os.path.join(OUTPUT_DIR, "top20_gene_priority_score.png")
plt.savefig(plot1_path, dpi=300)
plt.show()

print("📊 Grafik kaydedildi:", plot1_path)

# ------------------------------------------------------------

# Grafik 2: Skor dağılımı
plt.figure(figsize=(10, 5))
plt.hist(df_sorted["gene_priority_score"], bins=50)
plt.title("Gene Priority Score Distribution")
plt.xlabel("Score")
plt.ylabel("Number of genes")
plt.tight_layout()

# Kaydet
plot2_path = os.path.join(OUTPUT_DIR, "gene_priority_score_distribution.png")
plt.savefig(plot2_path, dpi=300)
plt.show()

print("📊 Grafik kaydedildi:", plot2_path)
