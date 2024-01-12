## requires as input a bed file with 4 columns (unlike a bed file this file should not be 0-based)
## object 'bed' is a data.frame with 4 columns, with columns 1-3 all having to be numeric
## currently designed for hg19
## chromosome IDs as numbers without 'chr'

#setwd("~/arch_pipeline/")

library(lifecycle, lib.loc=c("/gpfs/space/home/danat95/R/x86_64-pc-linux-gnu-library/4.1/"))
library(rlang, lib.loc=c("/gpfs/space/home/danat95/R/x86_64-pc-linux-gnu-library/4.1/"))
library(data.table)
library(dplyr)

## the functions and the loop is not written by me, so while it works as exepected I'm good

## helper function
fi=function(x){
  y=strsplit(x,";")[[1]]
  y=c(as.numeric(strsplit(y[1],":")[[1]][2]),as.numeric(strsplit(y[length(y)],":")[[1]][2]))
  y
}

## function to calculate genetic distances from genetic map files
fa=function(x,m){
  map=m
  co=fi(x)
  out=NA
  if((co[2]-co[1])==0) out=0
  if(sum(map[,1]>co[1] & map[,1]<co[2])>0 & (co[2]-co[1])>0){
    p=which(map[,1]>co[1] & map[,1]<co[2])
    xx=map[(min(p)-1):(max(p)+1),]
    
    if(nrow(xx)>3){
        out=diff(xx[c(2,(nrow(xx)-1)),3])+
            (((xx[2,1]-co[1])/diff(xx[1:2,1]))*diff(xx[1:2,3]))+
            (((co[2]-xx[nrow(xx)-1,1])/diff(xx[(nrow(xx)-1):nrow(xx),1]))*diff(xx[(nrow(xx)-1):nrow(xx),3]))
    }
    if(nrow(xx)==3){
        out=((co[2]-xx[2,1])/diff(xx[2:3,1])*diff(xx[2:3,3]))+
            ((xx[2,1]-co[1])/diff(xx[1:2,1])*diff(xx[1:2,3]))
    }
    
  }
  if(sum(map[,1]>co[1])>0 & sum(map[,1]<co[2])>0 & (co[2]-co[1])>0 & sum(map[,1]>co[1] & map[,1]<co[2])==0){
    p=which(map[,1]>co[1])[1]
    xx=map[(min(p)-1):(max(p)),]
    out=diff(co)/(xx[2,1]-xx[1,1])*(xx[2,3]-xx[1,3])
  }
  out
}

args <- commandArgs(trailingOnly = TRUE)
#args <- c("22","deCode.hg38.txt","AAMap.hg38.txt")
name <- args[1]
#name <- 10 # 1241/1217
#name <- 22 # ~500


dt <- fread(paste0("res/chr",name,".haplotypes.R2_0_8.tsv"))
dt$CHROM <- gsub("X","23",dt$CHROM)

bed <- dt[,c(3,4,5,1)] %>% data.frame()
bed$CHROM <- as.numeric(gsub("^chr","",bed$CHROM)) 
#str(bed)

#chr <- 22
## code to calculate 
bed1=c()
for(chr in sort(unique(bed[,1]))){
  #cat(chr,"\n")
  b.sub=cbind(bed[bed[,1]==chr,,drop=F],NA,NA)
  ## deCode map
  map=data.frame(fread(args[2])) # deCode first
  map <- map[map$CHROM == chr,c(2,3,4)]
  ## African American map
  map1=data.frame(fread(args[3]))
  map1 <- map1[map1$CHROM == chr,c(2,3,4)]

  for(d in 1:nrow(b.sub)){
      b.sub[d,5]=fa(paste(chr,":",b.sub[d,2],";",chr,":",b.sub[d,3],sep=""),map)/(b.sub[d,3]-b.sub[d,2])*1000000
      b.sub[d,6]=fa(paste(chr,":",b.sub[d,2],";",chr,":",b.sub[d,3],sep=""),map1)/(b.sub[d,3]-b.sub[d,2])*1000000
  }  
  
  t=465000
  b.sub[,5]=sapply(1:nrow(b.sub),function(o) 1-pgamma(diff(unlist(b.sub[o,2:3])),
                                                      shape=2,rate=1/(1/((b.sub[o,5]*1e-8*t)/25))))
  b.sub[,6]=sapply(1:nrow(b.sub),function(o) 1-pgamma(diff(unlist(b.sub[o,2:3])),
                                                      shape=2,rate=1/(1/((b.sub[o,6]*1e-8*t)/25))))
  
  b.sub=cbind(b.sub,apply(b.sub[,2:3,drop=F],1,diff))
  bed1=rbind(bed1,b.sub)
}

bed1=cbind(bed1,p.adjust(bed1[,5],"BH"),p.adjust(bed1[,6],"BH"))

colnames(bed1)=c("CHROM","start","end","chr_st_end","P_ILS_deCode","P_ILS_AAmap","haplotype_length","FDR_ILS_deCode","FDR_ILS_AAmap")

## counting ILS

f <- bed1[(bed1$FDR_ILS_deCode < 0.05 & !is.na(bed1$FDR_ILS_deCode)) |
              (bed1$FDR_ILS_AAmap < 0.05 & !is.na(bed1$FDR_ILS_AAmap)),]


#sum(is.na(f$FDR_ILS_AAmap) & is.na(f$FDR_ILS_deCode))
#f[is.na(f$FDR_ILS_AAmap) & is.na(f$FDR_ILS_deCode),]

#sum(is.na(f$FDR_ILS_AAmap)) # 4 # 4
#sum(is.na(f$FDR_ILS_deCode)) # 1566 # 2034
#sum(is.na(f$FDR_ILS_AAmap) & (f$FDR_ILS_deCode > 0.05 & !is.na(f$FDR_ILS_deCode))) # 1 # 1
#sum(is.na(f$FDR_ILS_deCode) & (f$FDR_ILS_AAmap > 0.05 & !is.na(f$FDR_ILS_AAmap))) #1232 # 1841

#print(c(sd(bed1$haplotype_length),mean(bed1$haplotype_length),median(bed1$haplotype_length)))
#print(c(sd(f$haplotype_length),mean(f$haplotype_length),median(f$haplotype_length)))
#bed1[is.na(bed1$FDR_ILS_AAmap) & is.na(bed1$FDR_ILS_deCode),]


dt$ILS[dt$chr_st_end %in% f$chr_st_end] <- "ILS_PASS"
dt$ILS[is.na(dt$ILS)] <- "ILS_FAIL"
dt$ILS[dt$ILS == "ILS_FAIL" & dt$N_aSNP >= 10] <- "ILS_FAIL_but_10"
#table(dt$ILS)

bed1$CHROM <- gsub("^","chr",bed1$CHROM)
bed1$CHROM <- gsub("23","X",bed1$CHROM)
dt$CHROM <- gsub("23","X",dt$CHROM)
write.table(bed1,paste0("res/chr",name,".ILS_raw_values.tsv"), quote = F, row.names = F, sep = "\t")
write.table(dt[,c(1:8,12)],paste0("res/chr",name,".haplotypes.R2_0_8.ILS.tsv"), quote = F, row.names = F, sep = "\t")

