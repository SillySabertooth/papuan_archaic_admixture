i=$1
smpl=$2

pth2=/gpfs/space/home/danat95/pr3_Pap_French/archaic_genoms/four_hg38


bcftools view -S $smpl \
  ~/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr$i.filtered.SNV_INDEL_SV_phased_panel.vcf.gz |\
  bcftools norm -m - \
  --fasta-ref /gpfs/space/databases/broadinstitute/references/hg38/v0/Homo_sapiens_assembly38.fasta \
  -Ob -o data/chr$i.splt.bcf && bcftools index data/chr$i.splt.bcf

### leave only snps  
bcftools view -v snps data/chr$i.splt.bcf -Ob -o data/chr$i.splt.snps.bcf && \
  bcftools index data/chr$i.splt.snps.bcf


bcftools merge $pth2/chr$i.4_genomes.hg38.fixed_ref.masked.bcf \
  data/chr$i.splt.snps.bcf |\
  bcftools norm -m + | bcftools norm -m - --fasta-ref /gpfs/space/databases/broadinstitute/references/hg38/v0/Homo_sapiens_assembly38.fasta \
  -Ob -o data/chr$i.Arch.out_target.tmp.bcf

bcftools view -v snps data/chr$i.Arch.out_target.tmp.bcf |\
  bcftools annotate -x INFO \
  -Ob -o data/chr$i.Arch.out_target.bcf

bcftools index data/chr$i.Arch.out_target.bcf