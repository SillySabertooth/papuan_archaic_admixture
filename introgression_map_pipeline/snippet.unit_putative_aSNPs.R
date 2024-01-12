#setwd("~/arch_pipeline/")

library(lifecycle, lib.loc=c("/gpfs/space/home/danat95/R/x86_64-pc-linux-gnu-library/4.1/"))
library(rlang, lib.loc=c("/gpfs/space/home/danat95/R/x86_64-pc-linux-gnu-library/4.1/"))
library(data.table)
library(dplyr)
library(tidyr)


#chr=22
chr <- commandArgs(trailingOnly = TRUE)

outgroup <- fread(paste0("res/chr",chr,".outgroup.info.tsv"), header = F)
colnames(outgroup) <- c("ID","AN_out","AC_out","YRI_status")


arch <- fread(paste0("res/chr",chr,".arch.info.tsv"), header = F)
arch_samples <- fread("Arch.samples", header = F)$V1
colnames(arch) <- c("ID","AN_arch","AC_arch",arch_samples,"dummy")
arch <- arch %>% select(-dummy)

target <- fread(paste0("res/chr",chr,".target.info.tsv"), header = F)
colnames(target) <- c("ID","AN_trgt","AC_trgt")

all <- merge(outgroup,target) %>% merge(.,arch)

f_remove <- all[all$AC_out != 0 & all$AN_out == all$AC_out & all$AN_trgt == all$AC_trgt,]

all_f <- all[!all$ID %in% f_remove$ID,]

all_f$AF_trgt <- all_f$AC_trgt/all_f$AN_trgt
all_f$AF_trgt_arch <- ifelse(all_f$YRI_status %in% c("1/1","1|1"),1-all_f$AF_trgt,all_f$AF_trgt)



# 
#sec_str <- all[all$AN_arch == all$AC_arch & all$AN_trgt == all$AC_trgt,] #chr22 77/5748 0,0134
# if use vcf and not gvcf, could be spurious cases where YRI ./. but all target is 1/1 and all called Arch 1/1. 
# Hope they will filter down by R2 to some extend 
# maybe I need to check it or remove/mark them explicitly;  #chr22 77/5748 0,0134




tmp_arch_slice <- all_f %>% select(all_of(arch_samples)) #%>% paste()
tmp_arch_slice <- lapply(tmp_arch_slice, function(x) sub("1/0", "0/1", x, fixed=TRUE)) %>% data.frame()


# homozyg attaching
# was significantly faster than having pain to use dynamic programming in r (for me)
homo_alt <- apply(tmp_arch_slice,1,function(x) as.numeric(
  table(x)[names(table(x)) == "1/1"]))
homo_alt[sapply(homo_alt, function(x) length(x)==0L)] <- 0
homo_alt <- unlist(homo_alt)

homo_ref <- apply(tmp_arch_slice,1,function(x) as.numeric(
  table(x)[names(table(x)) == "0/0"]))
homo_ref[sapply(homo_ref, function(x) length(x)==0L)] <- 0
homo_ref <- unlist(homo_ref)

all_f$arch_homo <- ifelse(all_f$YRI_status %in% c("1/1","1|1"),homo_ref,homo_alt)

# hetero
hetero <- apply(tmp_arch_slice,1,function(x) as.numeric(
  table(x)[names(table(x)) == "0/1"]))
hetero[sapply(hetero, function(x) length(x)==0L)] <- 0
all_f$arch_hetero <- unlist(hetero)

#all_f[,15:16] %>% table()

all_f <- separate(all_f,"ID",c("CHROM","POS","REF","ALT"), sep = ":", remove = F)

write.table(all_f,paste0("res/chr",chr,".putative_archaic_variants.tsv"), row.names = F, sep = "\t", quote = F)
write.table(all_f$ID,paste0("res/chr",chr,".putative_archaic_variants.list"), row.names = F, sep = "\t", quote = F, col.names = F)
