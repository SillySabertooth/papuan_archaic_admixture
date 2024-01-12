chrs=`echo $1 | sed "s/,/ /g"`
fst_chr=`echo $chrs | cut -f1 -d" "`

echo $chrs $fst_chr

mkdir -p res_concat/

for file in putative_archaic_variants.tsv archaic_variants.annotated_reg_info.tsv \
  aSNPs.R2_0_8.tsv haplotypes.R2_0_8.ILS.seq_sim.tsv \
  ILS_raw_values.tsv seq_sim.collapsed_annotated.tsv seq_sim.tsv \
  seq_sim.wide_debug.tsv shared_positions_R2_0.5.development.tsv; do
head -n1 res/chr$fst_chr.$file > res_concat/$file
for chrom in $chrs; do
tail -n+2 res/chr$chrom.$file >> res_concat/$file
done; done

for chrom in $chrs; do
cat res/chr$chrom.linked_haplo_R2_0.5.development.tsv >> res_concat/linked_haplo_R2_0.5.development.tsv
done

for file in putative_archaic_variants.tsv archaic_variants.annotated_reg_info.tsv \
  aSNPs.R2_0_8.tsv haplotypes.R2_0_8.ILS.seq_sim.tsv \
  ILS_raw_values.tsv seq_sim.collapsed_annotated.tsv seq_sim.tsv \
  seq_sim.wide_debug.tsv shared_positions_R2_0.5.development.tsv linked_haplo_R2_0.5.development.tsv; do
bgzip -c res_concat/$file > res_concat/$file.gz && rm res_concat/$file
done