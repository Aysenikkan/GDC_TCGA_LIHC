import os
import pandas as pd

# outputs klasörü yoksa oluştur
os.makedirs("outputs", exist_ok=True)

maf = pd.read_csv("merged_LIHC_MAF.csv", low_memory=False)

target_genes = [
    "TP53", "TERT", "CTNNB1",
    "ARID1A", "RB1", "AXIN1", "PTEN"
]

target_maf = maf[maf["Hugo_Symbol"].isin(target_genes)]

gene_counts = target_maf["Hugo_Symbol"].value_counts()

gene_counts.to_csv("outputs/target_gene_mutation_counts.csv")
