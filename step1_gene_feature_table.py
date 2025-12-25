# =========================================================
# STEP 1: Gene Feature Table Creation
# TCGA-LIHC somatic mutation data
# =========================================================

import pandas as pd
import os

# ---------------------------------------------------------
# 1) Çalışma dizinini ayarla (gerekirse)
# ---------------------------------------------------------
# os.chdir("D:/ALSU/GDC_TCGA_LIHC")

# ---------------------------------------------------------
# 2) Birleştirilmiş MAF dosyasını oku
# ---------------------------------------------------------
maf_path = "merged_LIHC_MAF.csv"
df = pd.read_csv(maf_path, low_memory=False)

print("MAF dosyası yüklendi.")
print("Toplam mutasyon sayısı (satır):", df.shape[0])
print("Toplam sütun sayısı:", df.shape[1])

# ---------------------------------------------------------
# 3) Gerekli sütunları seç
# ---------------------------------------------------------
required_cols = [
    "Hugo_Symbol",
    "Tumor_Sample_Barcode",
    "Variant_Classification",
    "IMPACT",
    "hotspot"
]

df = df[required_cols]

print("\nKullanılan sütunlar:")
print(df.columns.tolist())

# ---------------------------------------------------------
# 4) Toplam hasta sayısını hesapla
# ---------------------------------------------------------
total_patients = df["Tumor_Sample_Barcode"].nunique()
print("\nToplam hasta sayısı:", total_patients)

# ---------------------------------------------------------
# 5) Gen bazlı özet metrikleri hesapla
# ---------------------------------------------------------

# Toplam mutasyon sayısı (gen başına)
mutation_counts = df.groupby("Hugo_Symbol").size()

# Kaç farklı hastada mutasyon var
patient_counts = df.groupby("Hugo_Symbol")["Tumor_Sample_Barcode"].nunique()

# HIGH impact mutasyon sayısı
high_impact_counts = (
    df[df["IMPACT"] == "HIGH"]
    .groupby("Hugo_Symbol")
    .size()
)

# Hotspot mutasyon sayısı
hotspot_counts = (
    df[df["hotspot"] == True]
    .groupby("Hugo_Symbol")
    .size()
)

# ---------------------------------------------------------
# 6) Hepsini tek tabloda birleştir
# ---------------------------------------------------------
gene_features = pd.DataFrame({
    "n_mutations": mutation_counts,
    "n_patients": patient_counts,
    "n_high_impact": high_impact_counts,
    "hotspot_count": hotspot_counts
})

# NaN olanları 0 yap (örneğin hiç high-impact yoksa)
gene_features = gene_features.fillna(0)

# Oran hesapla
gene_features["high_impact_ratio"] = (
    gene_features["n_high_impact"] / gene_features["n_mutations"]
)

# Hasta frekansı (%)
gene_features["patient_frequency"] = (
    gene_features["n_patients"] / total_patients
)

# ---------------------------------------------------------
# 7) Sonuçları sırala (en çok mutasyona uğrayan genler üstte)
# ---------------------------------------------------------
gene_features = gene_features.sort_values(
    by="n_mutations",
    ascending=False
)

print("\nGen özet tablosu oluşturuldu.")
print("Toplam gen sayısı:", gene_features.shape[0])

print("\nİlk 10 gen:")
print(gene_features.head(10))

# ---------------------------------------------------------
# 8) Çıktıyı kaydet
# ---------------------------------------------------------
output_path = "outputs/gene_feature_table.csv"
gene_features.to_csv(output_path)

print("\nGen özet tablosu kaydedildi:")
print(output_path)
