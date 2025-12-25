import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# lifelines (survival analysis)
try:
    from lifelines import KaplanMeierFitter, CoxPHFitter
    from lifelines.statistics import logrank_test
except ImportError:
    raise ImportError(
        "lifelines yüklü değil. Kurmak için terminal/Anaconda Prompt:\n"
        "pip install lifelines"
    )

# ============================================================
# STEP 4B: Gene-based survival & recurrence analysis
# - Mutasyon (gene mutated vs not mutated) -> OS, DFS farkı var mı?
# Inputs:
#   outputs/clinical_prepared.csv
#   outputs/followup_prepared.csv
#   merged_LIHC_MAF.csv  (veya merged MAF dosyan)
# Optional:
#   outputs/gene_priority_score.csv (gene seçimini top N ile sınırlamak için)
# Outputs:
#   outputs/step4b_os_gene_results.csv
#   outputs/step4b_dfs_gene_results.csv
#   outputs/step4b_plots_os/*.png
#   outputs/step4b_plots_dfs/*.png
# ============================================================

BASE_DIR = r"D:\ALSU\GDC_TCGA_LIHC"
OUT_DIR = os.path.join(BASE_DIR, "outputs")
os.makedirs(OUT_DIR, exist_ok=True)

# ---- Input paths
CLIN_PATH = os.path.join(OUT_DIR, "clinical_prepared.csv")
FU_PATH   = os.path.join(OUT_DIR, "followup_prepared.csv")

# merged MAF yolu (senin dosyana göre güncelle)
MAF_PATH  = os.path.join(BASE_DIR, "merged_LIHC_MAF.csv")

# (opsiyonel) gen skor dosyası varsa top gen seçmek için
SCORE_PATH = os.path.join(OUT_DIR, "gene_priority_score.csv")

# ---- Params (istersen değiştir)
TOP_N_GENES = 500        # 14k genin hepsini yapmak ağır olabilir; önce 500 öneriyorum
MIN_MUT_PATIENTS = 10    # mutasyonlu grupta en az kaç hasta olsun
MIN_WT_PATIENTS  = 10    # mutasyonsuz grupta en az kaç hasta olsun
SAVE_TOP_PLOTS = 15      # en anlamlı kaç genin grafiğini kaydedelim (OS ve DFS ayrı)
ALPHA = 0.05

# ---- Output paths
OS_RES_PATH  = os.path.join(OUT_DIR, "step4b_os_gene_results.csv")
DFS_RES_PATH = os.path.join(OUT_DIR, "step4b_dfs_gene_results.csv")

PLOT_OS_DIR  = os.path.join(OUT_DIR, "step4b_plots_os")
PLOT_DFS_DIR = os.path.join(OUT_DIR, "step4b_plots_dfs")
os.makedirs(PLOT_OS_DIR, exist_ok=True)
os.makedirs(PLOT_DFS_DIR, exist_ok=True)

print("📥 Dosyalar okunuyor...")
clin = pd.read_csv(CLIN_PATH)
fu   = pd.read_csv(FU_PATH)
maf  = pd.read_csv(MAF_PATH, low_memory=False)

print("clinical_prepared:", clin.shape)
print("followup_prepared:", fu.shape)
print("merged MAF:", maf.shape)

# ------------------------------------------------------------
# 1) Patient ID eşleştirme
# clinical/followup 'patient_id' formatı: TCGA-XX-XXXX
# MAF Tumor_Sample_Barcode: TCGA-XX-XXXX-01A-... -> ilk 12 karakter patient
# ------------------------------------------------------------
required_maf_cols = ["Hugo_Symbol", "Tumor_Sample_Barcode"]
missing_maf = [c for c in required_maf_cols if c not in maf.columns]
if missing_maf:
    raise ValueError(f"MAF dosyasında eksik kolonlar: {missing_maf}")

maf["patient_id"] = maf["Tumor_Sample_Barcode"].astype(str).str.upper().str.slice(0, 12)
maf["Hugo_Symbol"] = maf["Hugo_Symbol"].astype(str)

# Klinik ID'leri normalize
clin["patient_id"] = clin["patient_id"].astype(str).str.upper().str.slice(0, 12)
fu["patient_id"]   = fu["patient_id"].astype(str).str.upper().str.slice(0, 12)

# Ortak hastalar
patients_os  = set(clin["patient_id"].unique())
patients_dfs = set(fu["patient_id"].unique())

print("\n👤 OS hastaları:", len(patients_os))
print("👤 DFS hastaları:", len(patients_dfs))

# ------------------------------------------------------------
# 2) Analiz edilecek gen listesini belirle
# - varsa gene_priority_score.csv içinden top N al
# - yoksa MAF'tan en çok hastada görülen top N al
# ------------------------------------------------------------
gene_list = None

if os.path.exists(SCORE_PATH):
    score_df = pd.read_csv(SCORE_PATH)
    if "Hugo_Symbol" in score_df.columns:
        gene_list = score_df["Hugo_Symbol"].astype(str).head(TOP_N_GENES).tolist()
        print(f"\n✅ Gen listesi gene_priority_score.csv içinden alındı: Top {TOP_N_GENES}")
    else:
        print("\n⚠ gene_priority_score.csv bulundu ama Hugo_Symbol yok, MAF'a düşüyorum...")

if gene_list is None:
    # MAF'tan gene patient count çıkar
    tmp = maf.groupby("Hugo_Symbol")["patient_id"].nunique().sort_values(ascending=False)
    gene_list = tmp.head(TOP_N_GENES).index.tolist()
    print(f"\n✅ Gen listesi MAF içinden seçildi: Top {TOP_N_GENES} (hasta sayısına göre)")

# ------------------------------------------------------------
# 3) Hasta->mutasyon matrisi için hızlı yapı
# gene -> set(patient_id)
# ------------------------------------------------------------
print("\n⚙ Gen->hasta setleri hazırlanıyor...")
gene_to_patients = maf.groupby("Hugo_Symbol")["patient_id"].apply(lambda s: set(s)).to_dict()

# ------------------------------------------------------------
# Yardımcı: KM plot kaydet
# ------------------------------------------------------------
def save_km_plot(time, event, mutated_mask, gene, out_png, title_prefix):
    kmf = KaplanMeierFitter()

    plt.figure(figsize=(8, 6))

    # Mutated
    kmf.fit(time[mutated_mask], event[mutated_mask], label=f"{gene} MUT (n={mutated_mask.sum()})")
    ax = kmf.plot()

    # WT
    kmf.fit(time[~mutated_mask], event[~mutated_mask], label=f"{gene} WT (n={(~mutated_mask).sum()})")
    kmf.plot(ax=ax)

    plt.title(f"{title_prefix}: {gene}")
    plt.xlabel("Gün")
    plt.ylabel("Sağkalım Olasılığı")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()

# ------------------------------------------------------------
# 4) OS Analizi (log-rank + Cox HR)
# ------------------------------------------------------------
print("\n🧬 OS analizi (gene mutated vs WT) başlıyor...")
os_results = []

os_df = clin.dropna(subset=["OS_time", "OS_event"]).copy()
os_df["OS_time"] = pd.to_numeric(os_df["OS_time"], errors="coerce")
os_df["OS_event"] = pd.to_numeric(os_df["OS_event"], errors="coerce")
os_df = os_df.dropna(subset=["OS_time", "OS_event"])

for gene in gene_list:
    mut_patients = gene_to_patients.get(gene, set())
    # Bu gene mutasyonu var mı?
    os_df["mut"] = os_df["patient_id"].apply(lambda pid: 1 if pid in mut_patients else 0).astype(int)

    n_mut = int(os_df["mut"].sum())
    n_wt = int((os_df["mut"] == 0).sum())

    if n_mut < MIN_MUT_PATIENTS or n_wt < MIN_WT_PATIENTS:
        continue

    t = os_df["OS_time"].values
    e = os_df["OS_event"].values
    m = os_df["mut"].values.astype(bool)

    # Log-rank test
    lr = logrank_test(t[m], t[~m], e[m], e[~m])
    p = float(lr.p_value)

    # Median survival (KM)
    kmf = KaplanMeierFitter()
    kmf.fit(t[m], e[m])
    med_mut = float(kmf.median_survival_time_) if kmf.median_survival_time_ is not None else np.nan
    kmf.fit(t[~m], e[~m])
    med_wt = float(kmf.median_survival_time_) if kmf.median_survival_time_ is not None else np.nan

    # Cox HR (tek değişken: mut)
    hr = np.nan
    try:
        cox_df = os_df[["OS_time", "OS_event", "mut"]].copy()
        cox_df.columns = ["T", "E", "mut"]
        cph = CoxPHFitter()
        cph.fit(cox_df, duration_col="T", event_col="E")
        hr = float(np.exp(cph.params_["mut"]))
    except Exception:
        hr = np.nan

    os_results.append({
        "gene": gene,
        "n_mut": n_mut,
        "n_wt": n_wt,
        "p_value": p,
        "cox_hr_mut_vs_wt": hr,
        "median_OS_mut_days": med_mut,
        "median_OS_wt_days": med_wt
    })

os_res = pd.DataFrame(os_results).sort_values("p_value").reset_index(drop=True)
os_res.to_csv(OS_RES_PATH, index=False)
print("✅ OS sonuçları kaydedildi:", OS_RES_PATH)
print("OS test edilen gen sayısı:", os_res.shape[0])
print("\nTop 10 (OS) en küçük p-value:")
print(os_res.head(10))

# OS plot kaydet (top)
print(f"\n🖼 OS için top {SAVE_TOP_PLOTS} KM grafiği kaydediliyor...")
for i, row in os_res.head(SAVE_TOP_PLOTS).iterrows():
    gene = row["gene"]
    mut_patients = gene_to_patients.get(gene, set())
    dfp = os_df.copy()
    dfp["mut"] = dfp["patient_id"].apply(lambda pid: 1 if pid in mut_patients else 0).astype(int)
    m = dfp["mut"].values.astype(bool)
    out_png = os.path.join(PLOT_OS_DIR, f"OS_KM_{i+1:02d}_{gene}.png")
    save_km_plot(dfp["OS_time"].values, dfp["OS_event"].values, m, gene, out_png, "Overall Survival (OS)")
print("✅ OS plotlar kaydedildi:", PLOT_OS_DIR)

# ------------------------------------------------------------
# 5) DFS/PFS Analizi (log-rank + Cox HR)
# ------------------------------------------------------------
print("\n🧬 DFS/PFS analizi (gene mutated vs WT) başlıyor...")
dfs_results = []

dfs_df = fu.dropna(subset=["DFS_time", "DFS_event"]).copy()
dfs_df["DFS_time"] = pd.to_numeric(dfs_df["DFS_time"], errors="coerce")
dfs_df["DFS_event"] = pd.to_numeric(dfs_df["DFS_event"], errors="coerce")
dfs_df = dfs_df.dropna(subset=["DFS_time", "DFS_event"])

for gene in gene_list:
    mut_patients = gene_to_patients.get(gene, set())
    dfs_df["mut"] = dfs_df["patient_id"].apply(lambda pid: 1 if pid in mut_patients else 0).astype(int)

    n_mut = int(dfs_df["mut"].sum())
    n_wt = int((dfs_df["mut"] == 0).sum())

    if n_mut < MIN_MUT_PATIENTS or n_wt < MIN_WT_PATIENTS:
        continue

    t = dfs_df["DFS_time"].values
    e = dfs_df["DFS_event"].values
    m = dfs_df["mut"].values.astype(bool)

    lr = logrank_test(t[m], t[~m], e[m], e[~m])
    p = float(lr.p_value)

    kmf = KaplanMeierFitter()
    kmf.fit(t[m], e[m])
    med_mut = float(kmf.median_survival_time_) if kmf.median_survival_time_ is not None else np.nan
    kmf.fit(t[~m], e[~m])
    med_wt = float(kmf.median_survival_time_) if kmf.median_survival_time_ is not None else np.nan

    hr = np.nan
    try:
        cox_df = dfs_df[["DFS_time", "DFS_event", "mut"]].copy()
        cox_df.columns = ["T", "E", "mut"]
        cph = CoxPHFitter()
        cph.fit(cox_df, duration_col="T", event_col="E")
        hr = float(np.exp(cph.params_["mut"]))
    except Exception:
        hr = np.nan

    dfs_results.append({
        "gene": gene,
        "n_mut": n_mut,
        "n_wt": n_wt,
        "p_value": p,
        "cox_hr_mut_vs_wt": hr,
        "median_DFS_mut_days": med_mut,
        "median_DFS_wt_days": med_wt
    })

dfs_res = pd.DataFrame(dfs_results).sort_values("p_value").reset_index(drop=True)
dfs_res.to_csv(DFS_RES_PATH, index=False)
print("✅ DFS/PFS sonuçları kaydedildi:", DFS_RES_PATH)
print("DFS test edilen gen sayısı:", dfs_res.shape[0])
print("\nTop 10 (DFS) en küçük p-value:")
print(dfs_res.head(10))

print(f"\n🖼 DFS için top {SAVE_TOP_PLOTS} KM grafiği kaydediliyor...")
for i, row in dfs_res.head(SAVE_TOP_PLOTS).iterrows():
    gene = row["gene"]
    mut_patients = gene_to_patients.get(gene, set())
    dfp = dfs_df.copy()
    dfp["mut"] = dfp["patient_id"].apply(lambda pid: 1 if pid in mut_patients else 0).astype(int)
    m = dfp["mut"].values.astype(bool)
    out_png = os.path.join(PLOT_DFS_DIR, f"DFS_KM_{i+1:02d}_{gene}.png")
    save_km_plot(dfp["DFS_time"].values, dfp["DFS_event"].values, m, gene, out_png, "Disease-Free / Progression-Free (DFS/PFS)")
print("✅ DFS plotlar kaydedildi:", PLOT_DFS_DIR)

# ------------------------------------------------------------
# 6) Mini özet
# ------------------------------------------------------------
print("\n====================")
print("STEP 4B BİTTİ ✅")
print("====================")
print("OS results :", OS_RES_PATH)
print("DFS results:", DFS_RES_PATH)
print("OS plots   :", PLOT_OS_DIR)
print("DFS plots  :", PLOT_DFS_DIR)

if os_res.shape[0] > 0:
    best = os_res.iloc[0]
    print(f"\n🏁 OS en anlamlı gen: {best['gene']} (p={best['p_value']:.3g}, HR={best['cox_hr_mut_vs_wt']})")

if dfs_res.shape[0] > 0:
    best = dfs_res.iloc[0]
    print(f"🏁 DFS en anlamlı gen: {best['gene']} (p={best['p_value']:.3g}, HR={best['cox_hr_mut_vs_wt']})")
