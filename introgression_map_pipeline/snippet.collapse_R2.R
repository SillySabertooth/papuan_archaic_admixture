#setwd("~/arch_pipeline/")

library(lifecycle, lib.loc=c("/gpfs/space/home/danat95/R/x86_64-pc-linux-gnu-library/4.1/"))
library(rlang, lib.loc=c("/gpfs/space/home/danat95/R/x86_64-pc-linux-gnu-library/4.1/"))
library(data.table)
library(dplyr)


###########################################################


# the function
#N <- 100
#list_haplo <- tmp_list
#r2check <- 0.4

create_haplot_table <- function(N, list_haplo, r2check){
  
  hap_tag <- names(list_haplo)[N]
  haplo <- unlist(list_haplo[N]) %>% as.numeric()
  
  # min R2 of united set
  
  min_R2 <- min(data$R2[data$POS1 %in% haplo & data$POS2 %in% haplo])
  
  # exploring SNPs around # not very good though, might be better to loon on any R2 up to some extend, or make it bigger
  # up
  snps_10_up <- posits[posits > max(haplo)][1:10]
  snps_10_up <- snps_10_up[!is.na(snps_10_up)]
  
  tmp_up <- data[data$POS1 %in% haplo & data$POS2 %in% snps_10_up & data$R2 >= r2check,]
  snps_10_up_PASS <- tmp_up$POS2 %>% unique() %>% length()
  
  # down
  
  snps_10_down <- rev(posits[posits < min(haplo)])[1:10]
  snps_10_down <- snps_10_down[!is.na(snps_10_down)]
  
  tmp_down <- data[data$POS1 %in% haplo & data$POS2 %in% snps_10_down & data$R2 >= r2check,]
  snps_10_down_PASS <- tmp_down$POS2 %>% unique() %>% length()
  
  # AF
  
  chunk_map <- mapNeand[mapNeand$POS %in% haplo,]
  af_median <- round(median(chunk_map$AF_trgt_arch),4)
  
  return(data.frame(haplo_tag = hap_tag,
                    chr_st_end = paste0(chr,".",min(haplo),".",max(haplo)),
                    N_aSNP=length(haplo),
                    CHROM=chr,
                    start = min(haplo),
                    end = max(haplo),
                    length = max(haplo) - min(haplo),
                    medianAF = af_median,
                    min_R2 = min_R2,
                    R2_10SNPs_nearby_DOWN = snps_10_down_PASS,
                    R2_10SNPs_nearby_UP = snps_10_up_PASS,
                    R2_10SNPs_nearby_TOTAL = snps_10_down_PASS+snps_10_up_PASS))
  
}


########################################################

name <- commandArgs(trailingOnly = TRUE)
#name <- 10 # 
#name <- 22 # ~500
chr <- paste0("chr",name)


# neand map
mapNeand <- fread(paste0("res/chr",name,".putative_archaic_variants.tsv")) %>% select(c(2:5,18))


data <- fread(paste0("aux_r2/chr",name,".hap.ld.gz"))
# if you need Dprime then do another dumb fix of the default vcftools bad file (data.table says that it is)
data <- data[,c(1,2,3,5)]
colnames(data) <- c("CHROM","POS1","POS2","R2")

data <- data[data$POS1 %in% mapNeand$POS & data$POS2 %in% mapNeand$POS,] #just to be on safe side

data <- rbind(data, data[,c(1,3,2,4)] %>% setnames(.,c("POS1","POS2"),c("POS2","POS1")))

posits <- c(data$POS1,data$POS2) %>% unique() %>% sort()



R2_ <- 0.8
# some time ago I run it for different thresholds
#for (R2_ in c(0.5,0.8)){

dt <- data[data$R2 >= R2_,]

dt$distance <- abs(dt$POS2-dt$POS1)
dt <- dt[dt$distance < 1e5,] # thumb rule treshold

tmp_posits <- c(dt$POS1,dt$POS2) %>% unique() %>% sort() %>% as.character()

dt$POS1 <- as.character(dt$POS1)
dt$POS2 <- as.character(dt$POS2)

tmp_list <- list()

tmp <- dt[dt$CHR == chr,]
# random was better to some extend, to choose a good haplotype, maybe remove even bigger than 1m

while (nrow(tmp) != 0){
  first_pos <- as.character(tmp_posits[1])
  tmp_list[first_pos] <- list(first_pos)
  
  while(nrow(tmp) != 0 & any(unlist(tmp_list[first_pos]) %in% tmp$POS1)){
    
    # attach to the list position that mathces position in the list
    tmp_list[first_pos] <- list(unique(c(unlist(tmp_list[first_pos]),
                                          tmp$POS2[tmp$POS1 %in% unlist(tmp_list[first_pos])])))
    
    # remove them from the table
    #tmp <- tmp[!tmp$POS1 %in% unlist(tmp_list[first_pos]),]
    tmp <- tmp[!tmp$POS2 %in% unlist(tmp_list[first_pos]),]
    tmp_posits <- tmp_posits[!tmp_posits %in% unlist(tmp_list[first_pos])]
    
  }
  
}



res <- lapply(1:length(tmp_list), function(x) create_haplot_table(x,tmp_list,(R2_ - 0.1))) %>% rbindlist()

## preparing output table for the collapsing step

R2_val <- gsub("\\.","_",R2_)

marker_table <- reshape2::melt(tmp_list)
colnames(marker_table) <- c("POS","haplo_tag")
#marker_table$haplo_tag <- gsub("^",paste0(name,"_"),marker_table$haplo_tag)
#marker_table$CHROM <- chr


marker_table <- merge(marker_table, res[,c(1,2,4)]) %>% select(chr_st_end,CHROM,POS)

write.table(res[,-1], paste0("res/chr",name,".haplotypes.R2_",R2_val,".tsv"), sep = "\t", row.names = F, quote = F)
write.table(marker_table, paste0("res/chr",name,".aSNPs.R2_",R2_val,".tsv"), sep = "\t", row.names = F, quote = F)



#let's say hotfix
#marker_table <- reshape2::melt(tmp_list)
#colnames(marker_table) <- c("aSNPs_hg38","haplo_tag_hg38")

haplo <- marker_table$chr_st_end[1]
linked_list <- list()
total_list <- list()

marker_table$POS <- as.numeric(marker_table$POS)

# stays here in case one will lauch it for 0.5 threshold
if (R2_ == 0.8){
  dt_0.5 <- data[data$R2 >= 0.5,]
  dt_0.5$distance <- abs(dt_0.5$POS2-dt_0.5$POS1)
  dt_0.5 <- dt_0.5[dt_0.5$distance < 1e5,] # thumb rule threshold here???
  links05 <- merge(dt_0.5, marker_table, by.x = c("CHROM","POS1"), by.y = c("CHROM","POS"))
  
  # get all linked to haplotypes
  
  linked_set <- merge(dt_0.5, marker_table, by.x = c("CHROM","POS1"), by.y = c("CHROM","POS")) %>% 
    merge(.,marker_table,by.x = c("CHROM","POS2"),by.y = c("CHROM","POS"))
  linked_set_tmp <- linked_set[linked_set$chr_st_end.x != linked_set$chr_st_end.y,]
  linked_set_tmp[,6:7] %>% unique() %>% nrow()
  
  write.table(unique(linked_set_tmp[,6:7]), paste0("res/",chr,".linked_haplo_R2_0.5.development.tsv"),col.names = F, sep = "\t", row.names = F, quote = F)
  
#  for (haplo in unique(marker_table$chr_st_end)){
  create_rare_shared <- function(haplo){
    set08 <- marker_table$POS[marker_table$chr_st_end == haplo]
    start <- min(set08)
    end <- max(set08)
    # positions linked to set aSNPs
    pos05 <- dt_0.5[dt_0.5$POS1 %in% set08 | dt_0.5$POS2 %in% set08,]
    pos05 <- unique(c(pos05$POS1,pos05$POS2))
    
    # positions linked to set08 but not present in it
    pos05 <- pos05[!pos05 %in% set08]
    pos05_other_sets <- links05[links05$POS2 %in% pos05,]
    pos05_other_sets$chr_st_end_target_set <- haplo
    
    pos05_other_sets <- pos05_other_sets[pos05_other_sets$POS2 >= start & pos05_other_sets$POS2 <= end,]
    pos05_other_sets <- pos05_other_sets[pos05_other_sets$chr_st_end != pos05_other_sets$chr_st_end_target_set,]
  }
  
  rareshared <- lapply(unique(marker_table$chr_st_end), function(x) create_rare_shared(x)) %>% rbindlist()
  rareshared[,6:7] %>% unique() %>% nrow()
  
  colnames(rareshared) <- c("CHROM","POS_belong_hap1","POS_linked_lay_notbelong_hap2_canbelong_hap1","R2","distance","hap1","hap2")
  
  write.table(unique(rareshared), paste0("res/chr",name,".shared_positions_R2_0.5.development.tsv"), sep = "\t", row.names = F, quote = F)
  
  
}


#}
