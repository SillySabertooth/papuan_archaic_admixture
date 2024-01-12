setwd("~/arch_pipeline_snakemake/")



# compare 1: maps

map_cur <- fread("../pr3_Pap_French/PAP_map.v5.yri_gvcf.arch_mask.tsv")


map_new <- fread("res_concat/putative_archaic_variants.tsv.gz")

sum(map_new$ID %in% map_cur$ID_PAP_hg38)

tst <- map_new[!map_new$ID %in% map_cur$ID_PAP_hg38,]
tst$AC_arch %>% table()
tst$AC_trgt[tst$AC_arch == 8] %>% table() # 3K, do they still here after collapsing?

clps <- fread("res_concat/aSNPs.R2_0_8.tsv.gz")
rst <- merge(clps,tst)
rst$AC_trgt[rst$AC_arch == 8] %>% table() # 3K, do they still here after collapsing?

#!# update map with regions and ils-pass,seq_sim
#!# add ID to some chr_pos?

# compare 2: regions

reg_cur <- fread("../pr3_Pap_French/index.ILS_counted.R2_0_8.map_upd.seq_sim.tsv")
reg_cur$chr_st_end <- gsub("^","chr",reg_cur$chr_st_end) 
#write.table(reg_cur[order(reg_cur$chr,reg_cur$start,reg_cur$end),c(2,4,5,13)],"compare_beds/reg_cur.bed",quote = F, row.names = F, col.names = F, sep = "\t")

reg_new <- fread("res_concat/haplotypes.R2_0_8.ILS.seq_sim.tsv.gz")
reg_new$CHROM <- gsub("^chr","",reg_new$CHROM)
reg_new$CHROM <- as.numeric(reg_new$CHROM)
#write.table(reg_new[order(reg_new$CHROM,reg_new$start,reg_new$end),c(3,4,5,1)],"compare_beds/reg_new.bed",quote = F, row.names = F, col.names = F, sep = "\t")

sum(reg_new$chr_st_end %in% reg_cur$chr_st_end)
sum(reg_cur$chr_st_end %in% reg_new$chr_st_end)

b <- reg_new[!reg_new$chr_st_end %in% reg_cur$chr_st_end,]
bb <- reg_cur[!reg_cur$chr_st_end %in% reg_new$chr_st_end,]

mrg <- merge(reg_cur,reg_new, by = "chr_st_end")

sum(mrg$N_aSNP.x == mrg$N_aSNP.y)

sum(mrg$length.y == mrg$length.x)

sum(mrg$ILS.y == mrg$ILS.x)

mrg %>% select(seq_sim,seq_sim_hap_filt) %>% table()
mrg %>% select(seq_sim_hap_filt) %>% table()
mrgx <- mrg[mrg$seq_sim_hap_filt == "seq_sim_FAIL" & mrg$seq_sim != "FAIL",]

# all is good now
#	chr1.102907303.103249483
# chr1.12214862.12523691
# chr1.12228475.12475664

# no really need to do it
#module load bedtools2/2.30.0; module load htslib
#cd compare_beds/
#bedtools intersect -a reg_cur.bed -b reg_new.bed -wo > reg_cur_new.bed
#regs <- fread("compare_beds/")

# compare 3: seq_sim

seq_cur <- fread("../pr3_Pap_French/seq_sim_res_annotated_with_Skov_res.tsv.gz")
seq_cur$sample_haplo <- gsub('.{5}$', '', seq_cur$sample_haplo)
seq_cur <- unique(seq_cur)

#seq_cur2 <- seq_cur[,c(1:2)] %>% unique()
#seq_cur2 <- separate(seq_cur2,seq_sim_chr_st_end,c("chr","st","end"),sep = "\\.", remove = F)
#seq_cur2$chr <- gsub("^chr","",seq_cur2$chr)
#seq_cur2$chr <- as.numeric(seq_cur2$chr)
#seq_cur2$an <- paste(seq_cur2$sample_haplo,seq_cur2$seq_sim_chr_st_end,sep = ":")

#write.table(seq_cur2[order(seq_cur2$chr,seq_cur2$st,seq_cur2$end),c(3,4,5,6)],"compare_beds/seq_cur.bed",quote = F, row.names = F, col.names = F, sep = "\t")

seq_new <- fread("res_concat/seq_sim.tsv.gz")
seq_new <- seq_new[seq_new$haplo_type == "archaic_target",]
seq_new$haplo <- gsub('.{5}$', '', seq_new$haplo)
seq_new <- unique(seq_new)

seqs <- merge(seq_cur, seq_new, by.x = c("sample_haplo","seq_sim_chr_st_end"), by.y = c("haplo","sub_hap"))

cor.test(seqs$seq_sim_AltaiNeandertal,seqs$AltaiNeandertal.seq_sim)
cor.test(seqs$seq_sim_Vindija33.19,seqs$Vindija33.19.seq_sim)
cor.test(seqs$`seq_sim_Chagyrskaya-Phalanx`,seqs$`Chagyrskaya-Phalanx.seq_sim`)
cor.test(seqs$seq_sim_Denisova,seqs$Denisova.seq_sim)

seqs_3 <- seqs[seqs$hap_full == "chr3.5639699.5669003",]
seq_new_3 <- seq_new[seq_new$hap_full == "chr3.5639699.5669003",]

df <- fread("../pr3_Pap_French/seq_sim/res.09.03.annotated.tsv.gz")
df <- df[df$haplo_type == "archaic",]
dfb <- df[df$chr_st_end == "3.5639699.5669003",]
dtb <- dt[dt$hap_full == "3.5639699.5669003",]
#df$chr_st_end %>% 

#ggplot(seqs, aes(seqs$seq_sim_AltaiNeandertal,seqs$AltaiNeandertal.seq_sim))+
#  geom_point()

seq_new$haplo_type %>% table()
sq_wide <- fread("res_concat/seq_sim.wide_debug.tsv.gz")
sq_wide <- sq_wide[sq_wide$haplo_type == "archaic_target",]
#module load bedtools2/2.30.0; module load htslib
#cd compare_beds/
#bedtools intersect -a seq_cur.bed -b seq_new.bed -wo > seq_cur_new.bed && bgzip -c seq_cur_new.bed > seq_cur_new.bed.gz && rm seq_cur_new.bed

#sbatch launhc_isec_filer.sh
#cd isec_res_17.03/
#cat chunks.??.splt.tsv >> $att_name.isec.Skov_mine.tsv && rm chunks.??.splt*
#cat $att_name.isec.Skov_mine.tsv | wc -l; cat $att_name.res | wc -l

#bgzip -c $att_name.isec.Skov_mine.tsv > $att_name.isec.Skov_mine.tsv.gz && rm $att_name.isec.Skov_mine.tsv
#bgzip -c $att_name.res > $att_name.isec.Skov_mine.all.tsv.gz && rm $att_name.res


#4 loook into combined
cmbs <- fread("res_concat/seq_sim.collapsed_annotated.tsv.gz")
#!#cmbs[cmbs == Inf] <- NA # yes

cmbs2 <- cmbs[cmbs$haplo_type == "archaic_target" & cmbs$ILS != "ILS_FAIL",] 
cmbs$hap_full <- gsub("^chr","",cmbs$hap_full)

cmbs_cur <- fread("../pr3_Pap_French/seq_sim/res.09.03.collapsed.annotated.arhaics_too.tsv.gz")
cmbs_cur <- cmbs_cur[cmbs_cur$haplo_type == "archaic",]

cmbs_mer <- merge(cmbs,cmbs_cur, by.x = "hap_full", by.y = "chr_st_end")

ggplot(cmbs_mer, aes(cmbs_mer$AltaiNeandertal.seq_sim_median,cmbs_mer$AltaiNeandertal))+
  geom_point()

cor.test(cmbs_mer$Denisova.seq_sim_median,cmbs_mer$Denisova)
