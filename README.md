# BLASTn Comparative Analysis: query-wheat vs reference-wheat High/Low-Confidence Genes

This workflow runs **BLASTn (megablast)** to compare high-confidence (HC)/low-confidence coding sequences (CDS) between wheat cultivar **query-wheat** and the **reference-wheat**. It identifies shared and unique HC genes between assemblies.

---

## 1. SLURM Job Script (Atlas HPC)

```bash
#!/bin/bash
#SBATCH -N 1
#SBATCH -n 48
#SBATCH -p atlas
#SBATCH --mem=376GB
#SBATCH --time=72:00:00
#SBATCH -J blast_query-reference-wheat
#SBATCH -A genolabswheatphg

ml blast-plus/2.14.1

cd /where/this/saved/Ruby/8-hc-lc-genes/

blastn -task megablast \
  -query /where/this/saved/8-hc-lc-genes/WHEAT-QUERY.high.cds.fa \
  -subject /90daydata/oat_renewal/Ruby/8-hc-lc-genes/WHEAT_REFERENCE_HC_cds.fasta \
  -evalue 1e-20 -perc_identity 90 -max_hsps 1 -num_threads 16 \
  -outfmt '6 qseqid sseqid pident length qlen slen evalue bitscore' \
  > wheat-query_to_wheat-reference-hc.megablast.tsv
```

## 2. Filtering Alignments
```bash
awk 'BEGIN{OFS="\t"}{
  pid=$3; L=$4; ql=$5; sl=$6; e=$7;
  covq=100*L/ql; covs=100*L/sl;
  if (pid>=90 && covq>=70 && covs>=70 && e<=1e-20) print $0, covq, covs
}' querywheat_to_subjectwheat-hc.megablast.tsv > querywheat_to_subjectwheat-hc.present.tsv
```

## 3. Extract Gene IDs
```bash
# All gene IDs
grep '^>' WHEAT-QUERY.high.cds.fa | sed 's/^>//' \
  > querywheat-hc.ids

grep '^>' WHEAT_REFERENCE_HC_cds.fasta | sed 's/^>//; s/ .*$//' \
  > subjectwheat-hc.ids

# Wheat query genes present in wheat subject
cut -f1 querywheat_to_subjectwheat-hc.present.tsv | sort -u > querywheat.present_in_subjectwheat-hc.ids

# Wheat subject genes present in query wheat
cut -f2 querywheat_to_subjectwheat-hc.present.tsv | sort -u > subjectwheat.present_in_querywheat-hc.ids

# Unique gene sets
comm -23 <(sort querywheat-hc.ids) <(sort querywheat.present_in_subjectwheat-hc.ids) > querywheat-hc.unique.ids
comm -23 <(sort subjectwheat-hc.ids)     <(sort subjectwheat.present_in_querywheat-hc.ids) > subjectwheat-hc.unique.ids
```

## 4. Summary Report
```bash
printf "querywheat total\t%s\n"  "$(wc -l < querywheat-hc.ids)"
printf "querywheat shared\t%s\n" "$(wc -l < querywheat.present_in_subjectwheat-hc.ids)"
printf "querywheat unique\t%s\n" "$(wc -l < querywheat-hc.unique.ids)"
printf "subjectwheat total\t%s\n"     "$(wc -l < subjectwheat-hc.ids)"
printf "subjectwheat shared\t%s\n"    "$(wc -l < subjectwheat.present_in_querywheat-hc.ids)"
printf "subjectwheat unique\t%s\n"    "$(wc -l < subjectwheat-hc.unique.ids)"
```

Maintainer:


Ruby Mijan





