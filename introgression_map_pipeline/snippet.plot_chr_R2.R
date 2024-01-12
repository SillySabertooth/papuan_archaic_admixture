#setwd("~/arch_pipeline/")

library(lifecycle, lib.loc=c("/gpfs/space/home/danat95/R/x86_64-pc-linux-gnu-library/4.1/"))
library(rlang, lib.loc=c("/gpfs/space/home/danat95/R/x86_64-pc-linux-gnu-library/4.1/"))
library(data.table)
library(dplyr)

library(ggplot2)
library(grid)
library(gridExtra)


name <- commandArgs(trailingOnly = TRUE)

#name <- 22

# preprocess the data, add remove one triangle
data <- fread(paste0("aux_r2/chr",name,".hap.ld.gz"))

data <- data[,c(1,2,3,5)]
colnames(data) <- c("CHROM","POS1","POS2","R2")

data <- rbind(data, data[,c(1,3,2,4)] %>% setnames(.,c("POS1","POS2"),c("POS2","POS1")))


data$id <- ifelse(data$POS1 < data$POS2, paste(data$POS1,data$POS2, sep = "_"),paste(data$POS2,data$POS1, sep = "_"))
data <- data[order(data$POS1),]



data <- data %>% 
  group_by(id) %>%
  distinct(id, .keep_all = T)
data <- data[order(data$POS1),]

# try to plot additional passed aSNPs
r2_collapsed <- fread(paste0("res/chr",name,".aSNPs.R2_0_8.tsv"))
index_haplo <- fread(paste0("res/chr",name,".haplotypes.R2_0_8.ILS.seq_sim.tsv"))
index_haplo <- index_haplo$chr_st_end[index_haplo$ILS != "ILS_FAIL" & !is.na(index_haplo$haplo_type)]
r2_collapsed <- r2_collapsed[r2_collapsed$chr_st_end %in% index_haplo,]

data_cond <- data[data$POS1 %in% r2_collapsed$POS,]

data <- rbind(data, data_cond[,c(1,3,2,4)] %>% setnames(.,c("POS1","POS2"),c("POS2","POS1")))

# the whole plot
#data$POS1 <- factor(data$POS1)
#data$POS2 <- factor(data$POS2)

#pdf(paste0("chr_plots/chr",name,".whole_chr.pdf")) # if you need put it in the same list, add p <- list(ggplot)
#ggplot(data = data, aes(x=factor(POS1), y=factor(POS2), fill=R2)) + 
#  geom_tile()+
#  ggtitle(paste0("chr",name,".whole_chr"))+
#  scale_fill_gradient('R2', limits=c(0, 1), breaks = c(0, 0.2, 0.5, 0.8, 1),  low = "yellow", high = "darkblue")+
#  theme(legend.position = "none", axis.text.x=element_blank(),
#        axis.text.y=element_blank(),
#        axis.title = element_blank(), axis.ticks = element_blank(),
#        plot.title = element_text(size = 8))
#dev.off()


####
# chunk plots

chr_pos <- c(data$POS1,data$POS2) %>% unique() %>% sort()

#length(chr_pos)%/%3000

chunk <- 3000

vec <- seq(1,length(chr_pos),chunk)
vec <- c(vec,length(chr_pos))


plots <- lapply(1:(length(vec)-1), function(i) {
  down_pos <- chr_pos[vec[i]]
  up_pos <- chr_pos[vec[i+1]]
  tmp <- data[data$POS1 >= down_pos & data$POS1 <= up_pos &
                data$POS2 >= down_pos & data$POS2 <= up_pos,]
  
  ggplot(data = tmp, aes(x=factor(POS1), y=factor(POS2), fill=R2)) + 
    geom_tile()+
    ggtitle(paste0("chr",name,"_chunk.",down_pos,"_",up_pos))+
    scale_fill_gradient('R2', limits=c(0, 1), breaks = c(0, 0.2, 0.5, 0.8, 1),  low = "yellow", high = "darkblue")+
    theme(legend.position = "none", axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.title = element_blank(), axis.ticks = element_blank(),
          plot.title = element_text(size = 8))}) 

ggsave(
  filename = paste0("chr_plots/chr",name,".chunks.pdf"), 
  plot = marrangeGrob(plots, nrow=1, ncol=1,top=NULL), 
  width = 10, height = 10
)
