import os
import gzip
import pandas as pd
from tqdm import tqdm

maf_dir = "maf_files"
all_maf = []

maf_files = [f for f in os.listdir(maf_dir) if f.endswith(".maf.gz")]

for maf in tqdm(maf_files):
    with gzip.open(os.path.join(maf_dir, maf), 'rt') as f:
        df = pd.read_csv(f, sep='\t', comment='#', low_memory=False)
        all_maf.append(df)

merged_maf = pd.concat(all_maf, ignore_index=True)

# analiz için dışa aktar
merged_maf.to_csv("merged_LIHC_MAF.csv", index=False)

print("MAF birleştirildi ve kaydedildi:")
print(merged_maf.shape)
