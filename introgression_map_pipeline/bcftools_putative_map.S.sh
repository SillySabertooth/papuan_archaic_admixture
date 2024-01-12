
file=$1
i=$2
arch=$3
outgroup=$4
target=$5
arch_AN_max=$6

#arch_AN_max=$((`cat $arch | wc -l` * 2)) # count max possible arch

### create a not fixed list to exclude

# if use vcf and not gvcf, could be spurious cases where YRI ./. but all target is 1/1 and all called Arch 1/1. Hope they will filter down by R2 to some extend 
# maybe I need to check it or remove/mark them explicitly;  #chr22 77/5748 0,0134

# in arhaics could be 1/0

bcftools view -S $outgroup $file |\
  bcftools filter -e 'AC == 0 | AC == AN' |\
  bcftools query -f'%ID\n' -o aux_lists/chr${i}.outgroup_not_fixed.IDs.list
  
### look for one where we have some difference, and try to work with ./., so exclude where we have missing

bcftools view --exclude ID==@aux_lists/chr${i}.outgroup_not_fixed.IDs.list \
  -S aux_lists/arch_outgroup.samples.tmp $file |\
  bcftools filter -e "AC == 0 | AC == AN && AN > $arch_AN_max" |\
  bcftools query -f'%ID\n' -o aux_lists/chr${i}.outgroup_fixed_diff_arch.IDs.list # by second command, remove 1/1 same fixed in archaics and outgroup
  
### look for one wherer we have diff with the targer

bcftools view --include ID==@aux_lists/chr${i}.outgroup_fixed_diff_arch.IDs.list \
  -S $target $file | bcftools filter -i "AC != 0" |\
  bcftools query -f'%ID\n' -o aux_lists/chr${i}.outgroup_fixed_diff_arch_and_target.not_full.list
  
  
  
# sublist for rare variants

bcftools view -S $outgroup $file |\
  bcftools filter -i 'AC == AN & AC != 0' |\
  bcftools query -f'%ID\n' -o aux_lists/chr${i}.out_fixed_alt.IDs.list
  
bcftools view --include ID==@aux_lists/chr${i}.out_fixed_alt.IDs.list \
  -S $arch $file | bcftools filter -i "AC != AN" |\
  bcftools query -f'%ID\n' -o aux_lists/chr${i}.out_fixed_alt_to_arch.IDs.list 

bcftools view --include ID==@aux_lists/chr${i}.out_fixed_alt_to_arch.IDs.list  \
  -S $target $file | bcftools filter -e "AC != 0" |\
  bcftools query -f'%ID\n' -o aux_lists/chr${i}.rare_fixed.list

cat aux_lists/chr${i}.outgroup_fixed_diff_arch_and_target.not_full.list \
  aux_lists/chr${i}.rare_fixed.list | sort -u > aux_lists/chr${i}.outgroup_fixed_diff_arch_and_target.list
  
#bcftools view -S res/arch_outgroup.samples.tmp data/chr${i}_phased.Arch.YRI.PAPs.snps.all.bcf | bcftools filter -i "AC != 0 && AC == AN && AN <= 8" | bcftools view -H | head

# preparing files to merge by R snippet

bcftools view --include ID==@aux_lists/chr${i}.outgroup_fixed_diff_arch_and_target.list \
  -S $outgroup $file |\
  bcftools query -f'%ID\t%AN\t%AC\t[%GT\t]\n' | cut -f1,2,3,4 > res/chr${i}.outgroup.info.tsv
  
bcftools view --include ID==@aux_lists/chr${i}.outgroup_fixed_diff_arch_and_target.list \
  -S $arch $file |\
  bcftools query -f'%ID\t%AN\t%AC\t[%GT\t]\n' -o res/chr${i}.arch.info.tsv
  

bcftools view --include ID==@aux_lists/chr${i}.outgroup_fixed_diff_arch_and_target.list \
  -S $target $file |\
  bcftools query -f'%ID\t%AN\t%AC\n' -o res/chr${i}.target.info.tsv