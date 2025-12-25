import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# STEP 3A: Skoru doğrulama + yorum raporu
# Girdi : outputs/gene_priority_score.csv
# Çıktı : outputs/step3A_top_genes.csv
#         outputs/step3A_known_driver_check.csv
#         outputs/step3A_report.txt
#         outputs/step3A_top20_score.png
#         outputs/step3A_score_vs_patientfreq.png
# ============================================================

BASE_DIR = r"D:\ALSU\GDC_TCGA_LIHC"   # <-- kendi yolun farklıysa değiştir
OUTPUT_DIR = os.path.join(BASE_DIR, "outputs")

INPUT_PATH = os.path.join(OUTPUT_DIR, "gene_priority_score.csv")

TOP_GENES_PATH = os.path.join(OUTPUT_DIR, "step3A_top_genes.csv")
DRIVER_CHECK_PATH = os.path.join(OUTPUT_DIR, "step3A_known_driver_check.csv")
REPORT_PATH = os.path.join(OUTPUT_DIR, "step3A_report.txt")

PLOT_TOP20_PATH = os.path.join(OUTPUT_DIR, "step3A_top20_score.png")
PLOT_SCATTER_PATH = os.path.join(OUTPUT_DIR, "step3A_score_vs_patientfreq.png")

os.makedirs(OUTPUT_DIR, exist_ok=True)

print("Okunan dosya:", INPUT_PATH)

# ------------------------------------------------------------
# 1) Skor tablosunu oku
# ------------------------------------------------------------
df = pd.read_csv(INPUT_PATH)

required = ["Hugo_Symbol", "gene_priority_score", "patient_frequency", "high_impact_ratio", "hotspot_ratio", "n_mutations", "n_patients"]
missing = [c for c in required if c not in df.columns]
if missing:
    raise ValueError(f"Eksik kolonlar var: {missing}\nMevcut kolonlar: {list(df.columns)}")

# numeric güvenliği
num_cols = ["gene_priority_score", "patient_frequency", "high_impact_ratio", "hotspot_ratio", "n_mutations", "n_patients"]
for c in num_cols:
    df[c] = pd.to_numeric(df[c], errors="coerce")
df = df.dropna(subset=["Hugo_Symbol"]).fillna(0)

# yeniden sırala (garanti)
df = df.sort_values("gene_priority_score", ascending=False).reset_index(drop=True)

print("Toplam gen:", df.shape[0])

# ------------------------------------------------------------
# 2) Top listeleri oluştur ve kaydet
# ------------------------------------------------------------
top50 = df.head(50).copy()
top100 = df.head(100).copy()

top_out = pd.concat(
    [top50.assign(top_k="top50"), top100.assign(top_k="top100")],
    ignore_index=True
)
top_out.to_csv(TOP_GENES_PATH, index=False)
print("✅ Top gen listeleri kaydedildi:", TOP_GENES_PATH)

# ------------------------------------------------------------
# 3) Basit 'bilinen driver' kontrolü (LIHC için sık geçenler)
# Not: Bu liste literatürde sık görülen genlerden "mini" bir kontrol listesidir.
# ------------------------------------------------------------
known_drivers = [
    "TP53", "CTNNB1", "AXIN1", "ARID1A", "ALB", "RB1",
    "TERT", "KEAP1", "NFE2L2", "APOB", "RPS6KA3",
    "ACVR2A", "BAP1", "CDKN2A", "PIK3CA"
]

df["is_known_driver_in_list"] = df["Hugo_Symbol"].isin(known_drivers).astype(int)

driver_hits = df[df["is_known_driver_in_list"] == 1].copy()
driver_hits["rank"] = driver_hits.index + 1
driver_hits = driver_hits[["Hugo_Symbol", "rank", "gene_priority_score", "patient_frequency", "high_impact_ratio", "hotspot_ratio", "n_mutations", "n_patients"]]
driver_hits = driver_hits.sort_values("rank")

driver_hits.to_csv(DRIVER_CHECK_PATH, index=False)
print("✅ Driver kontrol çıktısı:", DRIVER_CHECK_PATH)

# ------------------------------------------------------------
# 4) Grafikler (outputs'a kaydet)
# ------------------------------------------------------------

# Grafik 1: Top 20 barplot
top20 = df.head(20).copy()
plt.figure(figsize=(12, 6))
plt.bar(top20["Hugo_Symbol"], top20["gene_priority_score"])
plt.xticks(rotation=75, ha="right")
plt.title("STEP3A - Top 20 Gene Priority Score (LIHC)")
plt.xlabel("Gene")
plt.ylabel("Gene Priority Score")
plt.tight_layout()
plt.savefig(PLOT_TOP20_PATH, dpi=200)
plt.close()
print("✅ Grafik kaydedildi:", PLOT_TOP20_PATH)

# Grafik 2: Skor vs patient_frequency (scatter)
plt.figure(figsize=(8, 6))
plt.scatter(df["patient_frequency"], df["gene_priority_score"], s=10)
plt.title("STEP3A - Score vs Patient Frequency")
plt.xlabel("Patient Frequency")
plt.ylabel("Gene Priority Score")
plt.tight_layout()
plt.savefig(PLOT_SCATTER_PATH, dpi=200)
plt.close()
print("✅ Grafik kaydedildi:", PLOT_SCATTER_PATH)

# ------------------------------------------------------------
# 5) Kısa rapor yaz
# ------------------------------------------------------------
with open(REPORT_PATH, "w", encoding="utf-8") as f:
    f.write("STEP 3A REPORT - LIHC Gene Priority Score\n")
    f.write("========================================\n\n")
    f.write(f"Total genes: {df.shape[0]}\n")
    f.write(f"Top gene: {df.loc[0,'Hugo_Symbol']}  (score={df.loc[0,'gene_priority_score']:.4f})\n\n")

    f.write("Top 10 genes by score:\n")
    for i in range(10):
        g = df.loc[i, "Hugo_Symbol"]
        s = df.loc[i, "gene_priority_score"]
        pf = df.loc[i, "patient_frequency"]
        hi = df.loc[i, "high_impact_ratio"]
        hr = df.loc[i, "hotspot_ratio"]
        f.write(f"{i+1:02d}. {g:10s} score={s:.4f}  patient_freq={pf:.4f}  high_impact={hi:.4f}  hotspot={hr:.4f}\n")

    f.write("\nKnown driver mini-list hits (gene, rank):\n")
    if driver_hits.shape[0] == 0:
        f.write("No hits from the mini driver list.\n")
    else:
        for _, row in driver_hits.iterrows():
            f.write(f"- {row['Hugo_Symbol']}  rank={int(row['rank'])}  score={row['gene_priority_score']:.4f}\n")

    f.write("\nInterpretation notes:\n")
    f.write("- If known drivers appear in top ranks, score is biologically plausible.\n")
    f.write("- If very large genes (e.g., TTN) appear too high, there may be gene-length bias.\n")
    f.write("- Next: STEP 3B will cluster genes using mutation features (unsupervised ML).\n")

print("\n✅ STEP 3A tamamlandı.")
print("Rapor:", REPORT_PATH)
print("Top gen csv:", TOP_GENES_PATH)
print("Driver check:", DRIVER_CHECK_PATH)
