#setwd("~/arch_pipeline_snakemake/")
#setwd("~/pr3_Pap_French/KG1_SPOP/EUR/")

library(lifecycle, lib.loc=c("/gpfs/space/home/danat95/R/x86_64-pc-linux-gnu-library/4.1/"))
library(rlang, lib.loc=c("/gpfs/space/home/danat95/R/x86_64-pc-linux-gnu-library/4.1/"))
library(dplyr)#, lib.loc=c("/gpfs/space/home/danat95/R/x86_64-pc-linux-gnu-library/4.1/"))
library(tidyr)#, lib.loc=c("/gpfs/space/home/danat95/R/x86_64-pc-linux-gnu-library/4.1/"))
library(data.table)

# here will be exlucde arhaics missings to count which positions are total in a region
# and for the rest, if they are not present in modern human vcf and out of analyzed region
# we will consider them as most probably homozygote reference
# if not gvcf feed into the analysis, it would be the best possible solution I guess

args_ <- commandArgs(trailingOnly = TRUE)
#args_ <- c("22","Arch.samples","aux_lists/out_random_10.RS_out.samples.tmp","PAP.samples",
#           "/gpfs/space/home/danat95/pr3_Pap_French/archaic_genoms/four_hg38/chr22.4_genomes.hg38.fixed_ref.masked.bcf") # many ./. in YRI

#args_ <- c("X","Arch.samples","aux_lists/out_random_10.RS_out.samples.tmp","EUR.samples",
#           "/gpfs/space/home/danat95/pr3_Pap_French/archaic_genoms/four_hg38/chrX.4_genomes.hg38.fixed_ref.masked.bcf")

#Rscript ../pipe_v3_23_10_23/snippet.seq_sim.R X ../Arch.samples aux_lists/out_random_10.RS_out.samples.tmp ../EUR.samples /gpfs/space/home/danat95/pr3_Pap_French/archaic_genoms/four_hg38/chrX.4_genomes.hg38.fixed_ref.masked.bcf


name <- args_[1]

#groups
arch_S <- fread(args_[2], header = F)$V1
outgroup_S <- fread(args_[3], header = F)$V1
target_S <- fread(args_[4], header = F)$V1

#cols
bcf_cols <- c("ID","CHROM","POS","REF","ALT")

# maps
papmap <- fread(paste0("res/chr",name,".putative_archaic_variants.tsv")) %>% 
  select(all_of(c(bcf_cols,"YRI_status","arch_hetero","arch_homo")))

haplo_index <- fread(paste0("res/chr",name,".aSNPs.R2_0_8.tsv"))
chr_df_aSNPs <- merge(haplo_index,papmap, by = c("CHROM","POS"))


regions <- fread(paste0("res/chr",name,".haplotypes.R2_0_8.ILS.tsv"))

#i <- 131
#i <- 100
final_res_table <- data.frame()
wide_res_table <- data.frame()
for (i in 1:nrow(regions)){
#for (i in 1:20){

# the whole idea is using aSNPs, find the proper phase of arhaic haplotypes, and extract them into specific borders from the putative region
# and then count seq sim in the whole region

reg_bcf <- paste0("chr",name,":",regions$start[i],"-",regions$end[i])
reg_name <- regions$chr_st_end[i]
## create files
# bcftools should be attached
#system(paste0("bcftools view -S arch_10out_target.samples.tmp -r ",reg_bcf," data_id/chr",name,"_phased.Arch.YRI.PAPs.snps.all.IDs.bcf |\
#              bcftools query -f'%ID\t%CHROM\t%POS\t%REF\t%ALT\t[%GT\t]\n' > seq_sim/tsv_tmp/",reg_name,".arhc_pap.all.tsv"))

system("mkdir -p seq_sim/tsv_tmp; mkdir -p seq_sim/empty_reg/")

system(paste0("bcftools view -S ",args_[2]," -r ",reg_bcf," ",args_[5]," |\
              bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\t[%GT\t]\n' > seq_sim/tsv_tmp/",reg_name,".arch_all.tsv"))

system(paste0("bcftools view -S aux_lists/arch_10out_target.samples.tmp -r ",reg_bcf," data_id/chr",name,".original_vcf_with_IDs.bcf |\
              bcftools view -e 'AC == 0' | bcftools query -f'%ID\t%CHROM\t%POS\t%REF\t%ALT\t[%GT\t]\n' > seq_sim/tsv_tmp/",reg_name,".arhc_pap.vary.tsv"))


### attach the datasets
# all SNPs from the region

df <- fread(paste0("seq_sim/tsv_tmp/",reg_name,".arhc_pap.vary.tsv"))
colnames(df) <- c(bcf_cols,fread("aux_lists/arch_10out_target.samples.tmp", header = F)$V1,"dummy")

df_no_out <- df %>% select(all_of(c(bcf_cols,arch_S,target_S))) 


## 1. extract aSNPs, and found target and reference haplotypes in the given samples 

chr_df_aSNPs_reg <- chr_df_aSNPs[chr_df_aSNPs$chr_st_end == reg_name,]
df_aSNPs <- merge(chr_df_aSNPs_reg,df_no_out, by = bcf_cols)

# split on haplotypes, parse 1/1 YRI

df_aSNPs_adj <- df_aSNPs
# since we put in the begging our target to the end of df
# we can get first sample index by this command
n <- which(startsWith(colnames(df_aSNPs_adj),target_S[1]) == TRUE) 

# for chrX, do pseudohaploid calls
#if (name == "X"){
#  Xtarget <- df %>% select(all_of(target_S))
#  vec_Male <- lapply(Xtarget, function(x) ifelse("1" %in% x | "0" %in% x,T,F)) %>% unlist()
#  males <- colnames(Xtarget)[vec_Male]
 ### #Xtarget_male <- Xtarget %>% select(all_of(colnames(Xtarget)[vec_Male]))
  
#  df_aSNPs_adj[,n:length(df_aSNPs_adj)] <- lapply(df_aSNPs_adj[,n:length(df_aSNPs_adj)], function(x) as.character(x))
#  df_aSNPs_adj[,n:length(df_aSNPs_adj)] <- lapply(df_aSNPs_adj[,n:length(df_aSNPs_adj)], function(x) plyr::revalue(x, c("1" = "1|1","0" = "0|0"),warn_missing = F))
  
#}


for (trgt in target_S){
  df_aSNPs_adj <- df_aSNPs_adj %>% separate(trgt,
                into = c(paste0(trgt,"_hap1"),
                         paste0(trgt,"_hap2")),
                sep = "\\|")
}



yri_alt <- df_aSNPs_adj[df_aSNPs_adj$YRI %in% c("1|1","1/1"),]

yri_alt[,n:length(df_aSNPs_adj)] <- lapply(yri_alt[,n:length(df_aSNPs_adj)], factor)
yri_alt[,n:length(df_aSNPs_adj)] <- lapply(yri_alt[,n:length(df_aSNPs_adj)], function(x) plyr::revalue(x, c("1" = "0","0" = "1"),warn_missing = F))

df_aSNPs_adj <- rbind(df_aSNPs_adj[!df_aSNPs_adj$YRI %in% c("1|1","1/1"),],
                      yri_alt)
df_aSNPs_adj[,n:length(df_aSNPs_adj)] <- lapply(df_aSNPs_adj[,n:length(df_aSNPs_adj)], as.character)
df_aSNPs_adj[,n:length(df_aSNPs_adj)] <- lapply(df_aSNPs_adj[,n:length(df_aSNPs_adj)], as.numeric)

df_aSNPs_adj <- df_aSNPs_adj[order(df_aSNPs_adj$POS),]




# find arch and sample same number of mod haplotypes
hap_sum <- apply(df_aSNPs_adj[,n:length(df_aSNPs_adj)], 2, sum)
haplo_names <- colnames(df_aSNPs_adj)[n:length(df_aSNPs_adj)]

hap_arch <- haplo_names[hap_sum >= nrow(df_aSNPs)/2] #removing those that have less 50%?

if (length(hap_arch) == 0){
  write.table(hap_arch,paste0("seq_sim/empty_reg/",reg_name,".seq_sim.EMPTY.tsv"),row.names = F, sep = "\t", quote = F, col.names = F)
  next
}


hap_mod <- haplo_names[!(haplo_names %in% hap_arch)]

df_aSNPs_arch <- df_aSNPs_adj %>% select(all_of(hap_arch))

#hap_c <- 4
#hap_c <- 294
# double indecies, sweeet :()
count_borders <- function(hap_c){
  
  N_1 <- which(df_aSNPs_arch[[hap_c]] == 1)
  #N_1 <- c(1,2,5,6)
  med_N_1 <- round(median(N_1))
  
  
  if (any(diff(N_1) > 2)){
    first_down <- rev(N_1[N_1 < med_N_1])[1]
    first_down.INDEX <- which(N_1 == first_down)
    
    first_up <- N_1[N_1 > med_N_1][1]
    first_up.INDEX <- which(N_1 == first_up)
    
    if (first_up != max(N_1)){
    
    while (N_1[first_up.INDEX] - N_1[first_up.INDEX-1] < 3 & first_up.INDEX != length(N_1)){
      first_up.INDEX <- first_up.INDEX+1
    } 
      first_up <- N_1[first_up.INDEX-1] # stupied workaround to the counter, eh
      }
    
    if (first_down != min(N_1)){
    while (abs(N_1[first_down.INDEX] - N_1[first_down.INDEX+1]) < 3 & first_down.INDEX != 1){
      first_down.INDEX <- first_down.INDEX-1
    }
      if (first_down.INDEX != 1){
      first_down <- N_1[first_down.INDEX+1] ## stupied workaround to the counter, eh
      } else {
        first_down <- N_1[first_down.INDEX]
      }
    }
    
  } else {
    first_up <- max(N_1)
    first_down <- min(N_1)
  }
  
  return(data.frame(haplo=colnames(df_aSNPs_arch)[hap_c],
                   start=df_aSNPs_adj$POS[first_down],
                   end=df_aSNPs_adj$POS[first_up],
                   str_N=first_down,
                   end_N=first_up))

}

# create a table for indivduals arch haplotypes' borders, do the assiginig of them to random modern
# so then we can count seq sim in borders for inferred haplotypes archaic, and for similar length of modern human
haplo_st_end <- lapply(1:length(df_aSNPs_arch), function(x) count_borders(x)) %>%
  do.call(rbind, .)
haplo_st_end$aSNPs_num <- haplo_st_end$end_N-haplo_st_end$str_N+1

#workaround to remove those that pass somehow upper check from the hap_arch list 
haplo_st_end <- haplo_st_end[haplo_st_end$aSNPs_num > 1 & haplo_st_end$aSNPs_num >= nrow(df_aSNPs)/2,]
hap_arch <- haplo_st_end$haplo %>% unique()

if (length(hap_arch) == 0){
  write.table(hap_arch,paste0("seq_sim/empty_reg/",reg_name,".seq_sim.EMPTY.tsv"),row.names = F, sep = "\t", quote = F, col.names = F)
  next
}

# add modern for comparison

if (length(hap_arch) < length(hap_mod)){
  set.seed(404)
  hap_mod <- sample(hap_mod,length(hap_arch))
  random_st_end <- haplo_st_end
  random_st_end$haplo <- hap_mod
  random_st_end[,4:6] <- NA
} else {
  random_st_end <- sample_n(haplo_st_end,length(hap_mod))
  random_st_end$haplo <- hap_mod
  random_st_end[,4:6] <- NA
}


haplo_coord <- rbind(haplo_st_end,random_st_end)
haplo_coord$sub_hap <- paste(unique(df$CHROM),haplo_coord$start,haplo_coord$end, sep = ".") 
haplo_coord$sub_length <- haplo_coord$end - haplo_coord$start

## adding outgroup, just the whole region

# create outgroup sample

outgroup_S_haplo <- c(gsub("$","_hap1",outgroup_S),gsub("$","_hap2",outgroup_S))
if (length(hap_arch) < length(outgroup_S_haplo)){
  set.seed(40)
  outgroup_S_haplo <- sample(outgroup_S_haplo,length(hap_arch))
}

outgroup_S_haplo_tbl <- data.frame(haplo=outgroup_S_haplo,
                                   start=regions$start[i],
                                   end=regions$end[i],
                                   sub_hap=regions$chr_st_end[i],
                                   sub_length=regions$length[i])
haplo_coord <- bind_rows(haplo_coord,outgroup_S_haplo_tbl)

### counting similaity

## prepraing the whole avaivable from vcf, maybe here is some places where out and target different, so it can worth to look through
## so if out is 0/0 and target is ./., than they should have different number of simil

#df_G <- fread(paste0("seq_sim/tsv_tmp/",reg_name,".arhc_pap.all.tsv"))
#colnames(df_G) <- c(bcf_cols,fread("arch_10out_target.samples.tmp", header = F)$V1,"dummy")
###df_no_out <- df %>% select(all_of(c(bcf_cols,arch_S,target_S))) 

##df_G <- df_G %>% select(all_of(c(bcf_cols,arch_S,choosen_samples,outgroup_S)))
#df_G_miss <- df_G[apply(df_G, 1, function(r) any(r == "./.")),]
##total_dfg <- nrow(df_G)


#### preparing the whole dataset
# all names and data subset
choosen_samples <- gsub('.{5}$', '', c(hap_arch,hap_mod)) %>% unique()

df_taget <- df %>% select(all_of(c(bcf_cols,arch_S,choosen_samples,outgroup_S)))

# for chrX, do pseudohaploid
#if (name == "X"){
  
#  Xcontrol <- df %>% select(all_of(c(choosen_samples,outgroup_S)))
#  vec_Male_control <- lapply(Xcontrol, function(x) ifelse("1" %in% x | "0" %in% x,T,F)) %>% unlist()
#  males_control <- colnames(Xcontrol)[vec_Male_control]
  
#  N_targ_X <- length(c(bcf_cols,arch_S))+1
#  df_taget[,N_targ_X:length(df_taget)] <- lapply(df_taget[,N_targ_X:length(df_taget)], function(x) as.character(x))
#  df_taget[,N_targ_X:length(df_taget)] <- lapply(df_taget[,N_targ_X:length(df_taget)], function(x) plyr::revalue(x, c("1" = "1|1","0" = "0|0"),warn_missing = F))
  
#}

for (tgrt_smpl in c(choosen_samples,outgroup_S)){
  df_taget <- df_taget %>% separate(tgrt_smpl,
                                            into = c(paste0(tgrt_smpl,"_hap1"),
                                                     paste0(tgrt_smpl,"_hap2")),
                                            sep = "\\|")
}


arch_change <- df_taget %>% select(all_of(arch_S))
arch_change  <- data.frame(lapply(arch_change, factor))
arch_change <- data.frame(lapply(arch_change, function(x) plyr::revalue(x, c("1/1" = "1","0/1"="0.5","1/0"="0.5","0/0" = "0"))))
arch_change  <- data.frame(lapply(arch_change, as.character))
colnames(arch_change) <- arch_S

df_taget <- df_taget %>% select(!all_of(arch_S)) %>% replace(is.na(.), "./.")
df_taget[df_taget == "./."] <- 0
df_taget <- cbind(arch_change,df_taget)

#####

df_taget <- df_taget %>% select(all_of(c(bcf_cols,arch_S,hap_arch,hap_mod,outgroup_S_haplo)))

# the function for counting

count_sim <- function(smpl){
  if (smpl %in% hap_arch){
    hap_type <- "archaic_target"
  } else if (smpl %in% outgroup_S_haplo) {
    hap_type <- "outgroup"
  } else if (smpl %in% hap_mod) {
    hap_type <- "modern_target"
  }

  tmp_df <- df_taget[df_taget$POS >= haplo_coord$start[haplo_coord$haplo == smpl] & 
                       df_taget$POS <= haplo_coord$end[haplo_coord$haplo == smpl],]

  sub_hap_ <- haplo_coord$sub_hap[haplo_coord$haplo == smpl]
  
  # comparing itself
  #arch <- "AltaiNeandertal"
  #arch <- "Vindija33.19"
  #arch <- "Chagyrskayza-Phalanx"
  #arch <- "Denisova"
  #arch <- "Chagyrskayza-Phalanx"
  #smpl <- "O517_DA_B00IEW1_hap2"
  #smpl <- "O517_DA_C001KLN_hap1"
  #smpl <- "NA18881_hap2"
  res_tmp <- list()

  for (arch in arch_S){
    tmp_df_not_miss <- tmp_df[!tmp_df[[arch]] %in%  c("0.5","./."),]
    
    sum_not_mathc <- sum(tmp_df_not_miss[[arch]] != tmp_df_not_miss[[smpl]]) + 
      sum(tmp_df[[arch]] == "0.5")/2 # here we sum all that doesn't match, as well positions that hetero; 
    #they will be same for all comparison with this arhc, while add differ between diff arch
    
    res_tmp[arch] <- list(c(N_variants=nrow(tmp_df),
                            no_miss_no_het=nrow(tmp_df_not_miss),
                            het=sum(tmp_df[[arch]] == "0.5"),
                            not_match=sum_not_mathc))
  }
  return(cbind(data.frame(haplo=smpl,
                          haplo_type=hap_type,
                          sub_hap=sub_hap_),
               t(data.frame(unlist(res_tmp)))
  ))
}

seq_sim_res <- lapply(c(hap_arch,hap_mod,outgroup_S_haplo), function(x) count_sim(x)) %>%
  do.call(rbind, .)



# attach the last table; from the scratch add here chrNum.st.end

seq_sim_res <- merge(haplo_coord[,c(1,4:6,8)], seq_sim_res)


seq_sim_res <- cbind(hap_full = regions$chr_st_end[i],
                     seq_sim_res)
seq_sim_res$aSNPs_full_reg[seq_sim_res$haplo_type == "archaic_target"] <- nrow(df_aSNPs)


############

#and the last step to get idea how many actual missings here


df_arch <- fread(paste0("seq_sim/tsv_tmp/",reg_name,".arch_all.tsv"))
colnames(df_arch) <- c(bcf_cols[-1],fread("Arch.samples", header = F)$V1,"dummy")
df_arch <- df_arch %>% select(-dummy)

haplo_coord_uniq <- haplo_coord[,c(7,2,3)] %>% unique()

count_miss <- function(N){
  res_tmp <- list()
  #arch <- "Denisova"
  
  # comparing itself
  for (arch in arch_S){
    
    # count missing
    tmp_df <- df_arch[df_arch$POS >= haplo_coord_uniq$start[N] & 
                        df_arch$POS <= haplo_coord_uniq$end[N],]
    
    res_tmp[arch] <- list(c(avail_sites=sum(tmp_df[[arch]] != "./.")))
  }
  return(cbind(data.frame(sub_hap=haplo_coord_uniq$sub_hap[N]),
               t(data.frame(unlist(res_tmp)))
  ))
}

miss_arch <- lapply(1:nrow(haplo_coord_uniq), function(x) count_miss(x)) %>%
  do.call(rbind, .)


seq_sim_res2 <- merge(seq_sim_res,miss_arch)

for (arch_ind in arch_S){
  seq_sim_res2[[paste0(arch_ind,".seq_sim")]] <- seq_sim_res2[[paste0(arch_ind,".not_match")]]/seq_sim_res2[[paste0(arch_ind,".avail_sites")]]
  seq_sim_res2[[paste0(arch_ind,".miss_ratio")]] <- seq_sim_res2[[paste0(arch_ind,".avail_sites")]]/seq_sim_res2$sub_length
  seq_sim_res2[[paste0(arch_ind,".var_ratio")]] <- seq_sim_res2[[paste0(arch_ind,".N_variants")]]/seq_sim_res2[[paste0(arch_ind,".avail_sites")]]
}


seq_sim_res3 <- seq_sim_res2 %>% select(all_of(c("hap_full","haplo","sub_hap","sub_length","aSNPs_num","haplo_type",
                                               gsub("$",".seq_sim",arch_S))))


# adding archaics for the comparison

count_match <- function(sample){
  res_tmp <- list()
  #sample <- "Chagyrskaya-Phalanx"
  #arch <- "Denisova"
  
  # comparing itself
  
  for (arch in arch_S){
    
    clean_db <- df_arch[!(df_arch[[arch]] %in% c("1/0","0/1","./.") | df_arch[[sample]] %in% c("1/0","0/1","./.")),]
    not_mtch <- sum(clean_db[[arch]] != clean_db[[sample]])
    
    hetero_mism <- sum((df_arch[[arch]] %in% c("1/0","0/1") | df_arch[[sample]] %in% c("1/0","0/1")) &
                          df_arch[[arch]] != df_arch[[sample]])/2
    
    len_cl <- sum(!(df_arch[[arch]] == "./." | df_arch[[sample]] == "./."))
    #clean_db[clean_db[[arch]] != clean_db[[sample]],]
    
    res_tmp[arch] <- list(c(seq_sim=(not_mtch+hetero_mism)/len_cl))
  }
  return(cbind(data.frame(haplo=sample),
               t(data.frame(unlist(res_tmp)))
  ))
}



arch_res <- lapply(arch_S, function(x) count_match(x)) %>%
  do.call(rbind, .)

arch_res <- cbind(data.frame(hap_full=regions$chr_st_end[i],
                                   sub_hap=regions$chr_st_end[i],
                                   sub_length=regions$length[i],
                             haplo_type="archaic"),
                  arch_res)

seq_sim_res4 <- bind_rows(seq_sim_res3,arch_res)
seq_sim_res4 <- seq_sim_res4[order(seq_sim_res4$haplo_type),]
seq_sim_res2 <- seq_sim_res2[order(seq_sim_res2$haplo_type),]
#seq_sim_res4[seq_sim_res4 == Inf] <- NA

#if chrX, remove seconde pseudo male haplotype

if (name == "X"){
  #second_male <- c(gsub("$","_hap2",males),gsub("$","_hap2",males_control)) %>% unique()
  second_male <- fread("res/males.tsv",header = F)$V1
  second_male <- c(gsub("$","_hap2",second_male))
  
  seq_sim_res2 <- seq_sim_res2[!(seq_sim_res2$haplo %in% second_male),]
  seq_sim_res4 <- seq_sim_res4[!(seq_sim_res4$haplo %in% second_male),]
}



final_res_table <- rbind(final_res_table,seq_sim_res4)
wide_res_table <- rbind(wide_res_table,seq_sim_res2)
system(paste0("rm seq_sim/tsv_tmp/",reg_name,".arch_all.tsv seq_sim/tsv_tmp/",reg_name,".arhc_pap.vary.tsv"))
}

# save
write.table(final_res_table,paste0("res/chr",name,".seq_sim.tsv"),row.names = F, sep = "\t", quote = F)
write.table(wide_res_table,paste0("res/chr",name,".seq_sim.wide_debug.tsv"),row.names = F, sep = "\t", quote = F)

#if (name == "X"){
#write.table(males,paste0("res/males.tsv"),row.names = F, sep = "\t", quote = F, col.names = F)
#}
