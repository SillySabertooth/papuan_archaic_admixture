#setwd("~/arch_pipeline/")

library(lifecycle, lib.loc=c("/gpfs/space/home/danat95/R/x86_64-pc-linux-gnu-library/4.1/"))
library(rlang, lib.loc=c("/gpfs/space/home/danat95/R/x86_64-pc-linux-gnu-library/4.1/"))
library(dplyr)
library(tidyr)
library(data.table)

args_ <- commandArgs(trailingOnly = TRUE)
#args_ <- c("22","Arch.samples","PAP.samples","Denisova")
dt <- fread(paste0("res/chr",args_[1],".seq_sim.tsv"))

arch_S <- fread(args_[2], header = F)$V1

#wide_dbg <- fread(paste0("res/",args_[1],".seq_sim.wide_debug.tsv"))
#sum(wide_dbg$aSNPs_num == 1, na.rm = T)
#dt$aSNPs_full_reg[!is.na(dt$aSNPs_full_reg)] %>% table()

# add arch name
# add check for st==end just in case

#one simple adjustment
dt$haplo_type[dt$haplo_type == "archaic"] <- dt$haplo[dt$haplo_type == "archaic"]


#https://stackoverflow.com/questions/21644848/summarizing-multiple-columns-with-dplyr
#https://stackoverflow.com/questions/74647770/r-is-it-possible-to-use-median-and-na-rm-together-in-summarise-all

dt_b <- dt %>% group_by(hap_full,haplo_type) %>% summarise(N_haplos=n(),across(gsub("$",".seq_sim",arch_S), list(median = ~median(.,,na = TRUE))))
#dt_b <- dt_b %>% setnames(.,gsub("$",".seq_sim_1",arch_S),gsub("$",".seq_sim_median",arch_S))

table(dt_b$haplo_type)
#gsub("$",".seq_sim_median",arch_S)
#dt_b$closer_arch <- colnames(dt_b)[apply(dt_b,1,which.min)]

dt_b <- data.table(dt_b)
str(dt_b)

v <- apply(dt_b[dt_b$haplo_type == "archaic_target",] %>% 
             select(all_of(gsub("$",".seq_sim_median",arch_S))),1,function(x) which(x==min(x, na.rm = T)))
v_prob <- list()

#k <- 7
for (k in 1:length(v)){
  if (length(v[[k]]) == 1){
    v_prob[k] <- names(v[[k]])
  } else {
    if (paste0(args_[4],".seq_sim_median") %in% names(v[[k]])){
      v_prob[k] <- "Neand_Denisova"}
      else {
        v_prob[k] <- "Neand"
      }
  }
}

v_prob <- v_prob %>% unlist()
v_prob <- gsub(".seq_sim_median","_closer",v_prob)




dt_b$closer_arch[dt_b$haplo_type == "archaic_target"] <- v_prob
dt_b$closer_arch[dt_b$haplo_type != "archaic_target"] <- dt_b$haplo_type[dt_b$haplo_type != "archaic_target"]

dt_b$haplo_type %>% table()
dt_b$closer_arch %>% table()

dt_target <- dt_b[dt_b$haplo_type == "archaic_target",]

if (args_[1] == "X"){
  males <- fread("res/males.tsv", header = F)$V1
  AN_haplo <- (length(fread(args_[3], header = F)$V1)-length(males))*2+length(males)
  dt_target <- cbind(haplo_freq = dt_target$N_haplos/AN_haplo, 
                     dt_target) # ehehehe, sample list parced now :l
} else {
  dt_target <- cbind(haplo_freq = dt_target$N_haplos/(length(fread(args_[3], header = F)$V1)*2), 
                     dt_target) # ehehehe, sample list parced now :l
}





df1 <- fread(paste0("res/chr",args_[1],".haplotypes.R2_0_8.ILS.tsv"))

dt_res <- merge(df1,dt_target, by.x = "chr_st_end", by.y = "hap_full", all.x = T)

#dt_res$closer_arch %>% table()
#dt_res$closer_arch[dt_res$ILS != "ILS_FAIL"] %>% table()
#sum(dt_res$ILS != "ILS_FAIL")
dt_res$seq_sim_hap_filt[!is.na(dt_res$N_haplos)] <- "seq_sim_PASS"
dt_res$seq_sim_hap_filt[is.na(dt_res$N_haplos)] <- "seq_sim_FAIL"
#dt_res$seq_sim_hap_filt %>% table()

dt_b2 <- merge(dt_b, df1[,c(1,9)], by.y = "chr_st_end", by.x = "hap_full")
#dt_b2[,c(2,9)] %>% table()


snps <- fread(paste0("res/chr",args_[1],".aSNPs.R2_0_8.tsv"))
put_map <- fread(paste0("res/chr",args_[1],".putative_archaic_variants.tsv")) 
snps_put <- merge(snps,put_map, all.y = T) %>% merge(.,dt_res[,c(1,2,6,10,17,9,18)], by = "chr_st_end", all.x = T)


write.table(dt_res, paste0("res/chr",args_[1],".haplotypes.R2_0_8.ILS.seq_sim.tsv"), row.names = F, sep = "\t", quote = F)
write.table(snps_put, paste0("res/chr",args_[1],".archaic_variants.annotated_reg_info.tsv"), row.names = F, sep = "\t", quote = F)
write.table(dt_b2, paste0("res/chr",args_[1],".seq_sim.collapsed_annotated.tsv"), row.names = F, sep = "\t", quote = F)




#b <- fread("../pr3_Pap_French/seq_sim/res.09.03.collapsed.annotated.arhaics_too.tsv.gz")
#b <- b[startsWith(b$chr_st_end,"22."),]
#write.table(b,"test.umap.old.22.tsv", row.names = F, quote = F, sep = "\t")
#b <- b[b$haplo_type == "archaic"]
#b$haplo_type %>% table()
#b$prob %>% table()

#b <- b[b$haplo_type == "archaic" & b$ILS != "ILS_FAIL",]
#b$haplo_type %>% table()
#b$prob %>% table()
