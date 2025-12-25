import os
import pandas as pd
import numpy as np

# =========================
# STEP 4A (v2): Prepare clinical outcomes
# - OS (overall survival): time + event
# - DFS/PFS (recurrence/progression-free): time + event
#
# Key fix:
# - Use follow_ups.days_to_follow_up as censor time (very important!)
# =========================

BASE_DIR = r"D:\ALSU\GDC_TCGA_LIHC"
OUT_DIR = os.path.join(BASE_DIR, "outputs")
os.makedirs(OUT_DIR, exist_ok=True)

CLIN_PATH = os.path.join(BASE_DIR, "clinical.tsv")
FU_PATH   = os.path.join(BASE_DIR, "follow_up.tsv")

def pick_col(df, candidates):
    lower_map = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]
    return None

def find_like(df, keywords):
    cols = []
    for c in df.columns:
        lc = c.lower()
        ok = True
        for k in keywords:
            if k.lower() not in lc:
                ok = False
                break
        if ok:
            cols.append(c)
    return cols

# -------------------------
# 1) Read clinical.tsv
# -------------------------
clin = pd.read_csv(CLIN_PATH, sep="\t", low_memory=False)
print("\n[clinical.tsv] shape:", clin.shape)

id_col = pick_col(clin, [
    "cases.submitter_id",
    "demographic.submitter_id",
    "submitter_id",
    "case_submitter_id",
    "cases.case_id",
    "case_id",
    "patient_id"
])
if id_col is None:
    raise ValueError("clinical.tsv: patient id kolonu bulunamadı.")

vital_col = pick_col(clin, ["demographic.vital_status", "vital_status"])
dtd_col   = pick_col(clin, ["demographic.days_to_death", "days_to_death"])

if vital_col is None:
    raise ValueError("clinical.tsv: vital_status kolonu bulunamadı.")
if dtd_col is None:
    print("UYARI: days_to_death bulunamadı. OS event (death) zorlaşabilir.")

print("\nSeçilen kolonlar (clinical):")
print("id_col   :", id_col)
print("vital_col:", vital_col)
print("dtd_col  :", dtd_col)

clin["patient_id"] = clin[id_col].astype(str).str.upper()
clin[vital_col] = clin[vital_col].astype(str).str.lower()

# OS_event: dead=1 alive=0
clin["OS_event"] = clin[vital_col].map(lambda x: 1 if "dead" in x else 0)

if dtd_col is not None:
    clin[dtd_col] = pd.to_numeric(clin[dtd_col], errors="coerce")

# clinical tarafında OS_time (dead olanlar için days_to_death)
clin["OS_time_clin"] = np.nan
if dtd_col is not None:
    clin.loc[clin["OS_event"] == 1, "OS_time_clin"] = clin.loc[clin["OS_event"] == 1, dtd_col]

clin_small = clin[["patient_id", "OS_event", "OS_time_clin"]].drop_duplicates()

# -------------------------
# 2) Read follow_up.tsv
# -------------------------
fu = pd.read_csv(FU_PATH, sep="\t", low_memory=False)
print("\n[follow_up.tsv] shape:", fu.shape)

fu_id_col = pick_col(fu, [
    "cases.submitter_id",
    "submitter_id",
    "case_submitter_id",
    "cases.case_id",
    "case_id",
    "patient_id"
])
if fu_id_col is None:
    raise ValueError("follow_up.tsv: patient id kolonu bulunamadı.")

fu["patient_id"] = fu[fu_id_col].astype(str).str.upper()

# Censor time (takip süresi) -> follow_ups.days_to_follow_up
fu_follow_col = pick_col(fu, ["follow_ups.days_to_follow_up", "days_to_follow_up", "follow_up.days_to_follow_up"])
if fu_follow_col is None:
    raise ValueError("follow_up.tsv: days_to_follow_up bulunamadı. (Ama sende vardı, normalde bulmalı!)")

fu[fu_follow_col] = pd.to_numeric(fu[fu_follow_col], errors="coerce")

# Recurrence / progression time kolonları
rec_cols  = find_like(fu, ["days_to_recurrence"])
prog_cols = find_like(fu, ["days_to_progression"])

# Bazı dosyalarda duplicate kolon isimleri olabiliyor
rec_time_col  = rec_cols[0] if len(rec_cols) > 0 else None
prog_time_col = prog_cols[0] if len(prog_cols) > 0 else None

if rec_time_col is not None:
    fu[rec_time_col] = pd.to_numeric(fu[rec_time_col], errors="coerce")
if prog_time_col is not None:
    fu[prog_time_col] = pd.to_numeric(fu[prog_time_col], errors="coerce")

print("\nSeçilen kolonlar (follow_up):")
print("id_col          :", fu_id_col)
print("follow_up_time  :", fu_follow_col)
print("rec_time_col    :", rec_time_col)
print("prog_time_col   :", prog_time_col)

# -------------------------
# 3) Build OS_time using:
# - If DEAD: use clinical days_to_death
# - If ALIVE: use max follow_ups.days_to_follow_up
# -------------------------
fu_censor = fu.groupby("patient_id", as_index=False)[fu_follow_col].max()
fu_censor = fu_censor.rename(columns={fu_follow_col: "followup_time"})

os_df = clin_small.merge(fu_censor, on="patient_id", how="left")

os_df["OS_time"] = np.nan
# dead -> OS_time_clin
os_df.loc[os_df["OS_event"] == 1, "OS_time"] = os_df.loc[os_df["OS_event"] == 1, "OS_time_clin"]
# alive -> followup_time
os_df.loc[os_df["OS_event"] == 0, "OS_time"] = os_df.loc[os_df["OS_event"] == 0, "followup_time"]

os_df = os_df[["patient_id", "OS_time", "OS_event"]].dropna(subset=["OS_time"])
os_df["OS_time"] = pd.to_numeric(os_df["OS_time"], errors="coerce")
os_df = os_df.dropna(subset=["OS_time"]).drop_duplicates(subset=["patient_id"])

os_out = os.path.join(OUT_DIR, "clinical_prepared.csv")
os_df.to_csv(os_out, index=False)

print("\n✅ clinical_prepared.csv kaydedildi:", os_out)
print("   patients:", os_df["patient_id"].nunique())
print("   deaths  :", int(os_df["OS_event"].sum()))
print("   alive   :", int((os_df["OS_event"] == 0).sum()))

# -------------------------
# 4) Build DFS/PFS:
# event=1 if recurrence/progression exists
# else event=0 and time=followup_time
# For each patient:
# - event time = min(recurrence/progression)
# - censor time = max(days_to_follow_up)
# -------------------------
fu_work = fu[["patient_id", fu_follow_col]].copy()
fu_work = fu_work.rename(columns={fu_follow_col: "followup_time"})

fu_work["event_time"] = np.nan

if rec_time_col is not None:
    fu_work["rec_time"] = fu[rec_time_col]
else:
    fu_work["rec_time"] = np.nan

if prog_time_col is not None:
    fu_work["prog_time"] = fu[prog_time_col]
else:
    fu_work["prog_time"] = np.nan

# event_time: min of rec/prog if exists
fu_work["event_time"] = fu_work[["rec_time", "prog_time"]].min(axis=1, skipna=True)

# per patient: earliest event time
event_per_patient = fu_work.dropna(subset=["event_time"]).groupby("patient_id", as_index=False)["event_time"].min()
event_per_patient["DFS_event"] = 1
event_per_patient = event_per_patient.rename(columns={"event_time": "DFS_time"})

# per patient: censor time = max followup
censor_per_patient = fu_work.groupby("patient_id", as_index=False)["followup_time"].max()
censor_per_patient["DFS_event"] = 0
censor_per_patient = censor_per_patient.rename(columns={"followup_time": "DFS_time"})

dfs = pd.concat([event_per_patient, censor_per_patient], ignore_index=True)

# if patient has an event, keep event row
dfs = dfs.sort_values(["patient_id", "DFS_event"], ascending=[True, False])
dfs = dfs.drop_duplicates(subset=["patient_id"], keep="first")
dfs = dfs.dropna(subset=["DFS_time"])
dfs["DFS_time"] = pd.to_numeric(dfs["DFS_time"], errors="coerce")
dfs = dfs.dropna(subset=["DFS_time"])

dfs_out = os.path.join(OUT_DIR, "followup_prepared.csv")
dfs.to_csv(dfs_out, index=False)

print("\n✅ followup_prepared.csv kaydedildi:", dfs_out)
print("   patients:", dfs["patient_id"].nunique())
print("   events (rec/prog):", int(dfs["DFS_event"].sum()))
print("   censored:", int((dfs["DFS_event"] == 0).sum()))
