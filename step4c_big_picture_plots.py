import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# =========================
# STEP 4C: Big picture plots (robust)
# Inputs:
#  - outputs/step4b_os_gene_results.csv
#  - outputs/step4b_dfs_gene_results.csv
# Outputs:
#  - outputs/step4c_big_picture/*.png
# =========================

BASE_DIR = r"D:\ALSU\GDC_TCGA_LIHC"
OUT_DIR = os.path.join(BASE_DIR, "outputs")
PLOT_DIR = os.path.join(OUT_DIR, "step4c_big_picture")
os.makedirs(PLOT_DIR, exist_ok=True)

OS_PATH  = os.path.join(OUT_DIR, "step4b_os_gene_results.csv")
DFS_PATH = os.path.join(OUT_DIR, "step4b_dfs_gene_results.csv")

os_df = pd.read_csv(OS_PATH)
dfs_df = pd.read_csv(DFS_PATH)

print("OS columns:", list(os_df.columns))
print("DFS columns:", list(dfs_df.columns))

def pick_first_existing(df, candidates):
    """Return first candidate column that exists in df (case-insensitive)."""
    lower_map = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]
    return None

# Gene kolonu (genelde 'gene' olur)
gene_col_os  = pick_first_existing(os_df,  ["gene", "Hugo_Symbol", "hugo_symbol"])
gene_col_dfs = pick_first_existing(dfs_df, ["gene", "Hugo_Symbol", "hugo_symbol"])
if gene_col_os is None or gene_col_dfs is None:
    raise ValueError("Gene kolonu bulunamadı. OS/DFS CSV içindeki gen kolonunu kontrol et.")

# p-value kolonu (seninkinde farklı isim olabilir)
p_candidates = ["logrank_p", "p_value", "pvalue", "p", "pval", "logrank_pvalue", "logrank_p_value"]
p_col_os  = pick_first_existing(os_df, p_candidates)
p_col_dfs = pick_first_existing(dfs_df, p_candidates)
if p_col_os is None or p_col_dfs is None:
    raise ValueError(f"p-value kolonu bulunamadı. Aranan adaylar: {p_candidates}")

# HR kolonu (cox model sonucu)
hr_candidates = ["cox_hr_mut_vs_wt", "hr", "hazard_ratio", "cox_hr", "cox_hr_mutant_vs_wt"]
hr_col_os  = pick_first_existing(os_df, hr_candidates)
hr_col_dfs = pick_first_existing(dfs_df, hr_candidates)
if hr_col_os is None or hr_col_dfs is None:
    raise ValueError(f"HR kolonu bulunamadı. Aranan adaylar: {hr_candidates}")

print("\nSeçilen kolonlar:")
print("OS  gene:", gene_col_os, " p:", p_col_os, " HR:", hr_col_os)
print("DFS gene:", gene_col_dfs," p:", p_col_dfs," HR:", hr_col_dfs)

# gerekli kolonları normalize edip ortak isim verelim
os_df = os_df.rename(columns={gene_col_os: "gene", p_col_os: "p_value", hr_col_os: "HR"})
dfs_df = dfs_df.rename(columns={gene_col_dfs: "gene", p_col_dfs: "p_value", hr_col_dfs: "HR"})

# numeric dönüşüm
os_df["p_value"] = pd.to_numeric(os_df["p_value"], errors="coerce")
dfs_df["p_value"] = pd.to_numeric(dfs_df["p_value"], errors="coerce")
os_df["HR"] = pd.to_numeric(os_df["HR"], errors="coerce")
dfs_df["HR"] = pd.to_numeric(dfs_df["HR"], errors="coerce")

# yardımcı fonksiyonlar
def neglog10p(p):
    p = np.clip(p, 1e-300, 1.0)
    return -np.log10(p)

def log2hr(hr):
    hr = np.clip(hr, 1e-9, None)
    return np.log2(hr)

os_df["neglog10_p"] = neglog10p(os_df["p_value"])
dfs_df["neglog10_p"] = neglog10p(dfs_df["p_value"])
os_df["log2HR"] = log2hr(os_df["HR"])
dfs_df["log2HR"] = log2hr(dfs_df["HR"])


# ------------------------------------------------------------
# 1) Volcano: log2(HR) vs -log10(p)
# ------------------------------------------------------------
plt.figure(figsize=(10, 6))
plt.scatter(os_df["log2HR"], os_df["neglog10_p"])
plt.axvline(0, linestyle="--")
plt.axhline(-np.log10(0.05), linestyle="--")
plt.title("OS: Etki (log2(HR)) vs Anlamlılık (-log10 p)")
plt.xlabel("log2(HR) (mutant vs WT)")
plt.ylabel("-log10(p)")
plt.tight_layout()
plt.savefig(os.path.join(PLOT_DIR, "bigpic_os_volcano.png"), dpi=200)
plt.close()

plt.figure(figsize=(10, 6))
plt.scatter(dfs_df["log2HR"], dfs_df["neglog10_p"])
plt.axvline(0, linestyle="--")
plt.axhline(-np.log10(0.05), linestyle="--")
plt.title("DFS/PFS: Etki (log2(HR)) vs Anlamlılık (-log10 p)")
plt.xlabel("log2(HR) (mutant vs WT)")
plt.ylabel("-log10(p)")
plt.tight_layout()
plt.savefig(os.path.join(PLOT_DIR, "bigpic_dfs_volcano.png"), dpi=200)
plt.close()


# ------------------------------------------------------------
# 2) Top 10 bar: -log10(p)
# ------------------------------------------------------------
def top_bar(df, title, filename, topn=10):
    tmp = df.sort_values("p_value").head(topn).copy()
    plt.figure(figsize=(10, 5))
    plt.bar(tmp["gene"], neglog10p(tmp["p_value"]))
    plt.xticks(rotation=45, ha="right")
    plt.title(title)
    plt.xlabel("Gen")
    plt.ylabel("-log10(p)")
    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_DIR, filename), dpi=200)
    plt.close()

top_bar(os_df,  "OS: En anlamlı Top 10 gen (-log10 p)",  "bigpic_os_top10_p.png",  10)
top_bar(dfs_df, "DFS: En anlamlı Top 10 gen (-log10 p)", "bigpic_dfs_top10_p.png", 10)


# ------------------------------------------------------------
# 3) OS vs DFS log2(HR) scatter (ortak genler)
# ------------------------------------------------------------
merged = pd.merge(
    os_df[["gene", "HR", "p_value"]],
    dfs_df[["gene", "HR", "p_value"]],
    on="gene",
    suffixes=("_OS", "_DFS")
)

if merged.shape[0] > 0:
    merged["log2HR_OS"] = log2hr(merged["HR_OS"])
    merged["log2HR_DFS"] = log2hr(merged["HR_DFS"])

    plt.figure(figsize=(7, 7))
    plt.scatter(merged["log2HR_OS"], merged["log2HR_DFS"])
    plt.axvline(0, linestyle="--")
    plt.axhline(0, linestyle="--")
    plt.title("OS vs DFS: log2(HR) karşılaştırma (ortak genler)")
    plt.xlabel("log2(HR) - OS")
    plt.ylabel("log2(HR) - DFS")
    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_DIR, "bigpic_os_vs_dfs_log2hr_scatter.png"), dpi=200)
    plt.close()
else:
    print("Uyarı: OS ve DFS tablolarında ortak gen bulunamadı, scatter çizilmedi.")


# ------------------------------------------------------------
# 4) Yön matrisi (Top 10 OS + Top 10 DFS genleri)
# ------------------------------------------------------------
top_os = os_df.sort_values("p_value").head(10)["gene"].tolist()
top_dfs = dfs_df.sort_values("p_value").head(10)["gene"].tolist()
genes = list(dict.fromkeys(top_os + top_dfs))

mat = []
for g in genes:
    row = []
    r_os = os_df.loc[os_df["gene"] == g, "HR"]
    row.append(r_os.values[0] if len(r_os) else np.nan)
    r_dfs = dfs_df.loc[dfs_df["gene"] == g, "HR"]
    row.append(r_dfs.values[0] if len(r_dfs) else np.nan)
    mat.append(row)

mat = np.array(mat, dtype=float)
mat_log2 = log2hr(mat)

plt.figure(figsize=(6, max(4, 0.35 * len(genes))))
plt.imshow(mat_log2, aspect="auto")
plt.yticks(range(len(genes)), genes)
plt.xticks([0, 1], ["OS log2(HR)", "DFS log2(HR)"])
plt.title("Top genler: Mutasyon etkisi yön özeti (log2 HR)")
plt.colorbar(label="log2(HR)")
plt.tight_layout()
plt.savefig(os.path.join(PLOT_DIR, "bigpic_direction_matrix_log2hr.png"), dpi=200)
plt.close()

print("\n✅ STEP 4C bitti. Grafikler kaydedildi:")
print("->", PLOT_DIR)
