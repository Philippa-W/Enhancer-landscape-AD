### Imperial CHAS ####
library(CHAS)
library(RColorBrewer)
library(ggplot2)

setwd("~/OneDrive - University of Exeter/ATAC/paper/CHAS/rerun/FC/OLIG2")
setwd("~/Library/CloudStorage/OneDrive-UniversityofExeter/paper/CHAS/rerun/FC/OLIG2")
df <- read.csv("All_samples_counts_case copy.txt",sep = "\t", header = TRUE, skip = 1)
colnames(df) <- gsub(".mLb.clN.sorted.bam", "", colnames(df))
df2 <- df
df <- read.csv("All_samples_counts_cont copy.txt",sep = "\t", header = TRUE, skip = 1)
colnames(df) <- gsub("X.rds.general.user.pmwells.projects.epinott.live.user_analysed_data.Philippa.batch3.olig2.qc.cont.","", colnames(df))
colnames(df) <- gsub(".mLb.clN.sorted.bam", "", colnames(df))
df3 <- cbind(df2,df[-c(1:6)])
OLIG2 <- df3
write.csv(df3, "All_samples_counts_OLIG2.csv", row.names = FALSE)
setwd("~/Library/CloudStorage/OneDrive-UniversityofExeter/paper/CHAS/rerun/FC/PU1")
df <- read.csv("All_samples_counts_case copy.txt",sep = "\t", header = TRUE, skip = 1)
df <- read.csv("counts_qcd.csv")
PU1 <- df[-1]
setwd("~/Library/CloudStorage/OneDrive-UniversityofExeter/paper/CHAS")
NEUN <- read.csv("All_samples_counts_incsexchr.txt",sep = "\t", header = TRUE, skip = 1)
colnames(NEUN) <- gsub(".mLb.clN.sorted.bam","", colnames(NEUN))
setwd("~/Library/CloudStorage/OneDrive-UniversityofExeter/paper/CHAS/rerun")
write.csv(NEUN,"NEUN_fc.csv", row.names = FALSE)
write.csv(OLIG2,"OLIG2_fc.csv", row.names = FALSE)
write.csv(PU1,"PU1_fc.csv", row.names = FALSE)


cells <- c("PU1", "OLIG2", "NEUN")

dfmm <-data.frame(Metadata_per_file)
dfmm$AD <- ifelse(dfmm$Braak...Stage<3,0,1)
for (cell in 1:length(cells)){
  setwd("~/Library/CloudStorage/OneDrive-UniversityofExeter/paper/CHAS/rerun")
  peaksB <- read.csv(paste0(cells[cell],"_fc.csv"))
  peaksA <- peaksB[c(2,3,4,1)]
  row.names(peaksB) <- peaksB$Geneid
  peaksB <- peaksB[-c(1:6)]
  a <- as.data.frame(colnames(peaksB))
  colnames(peaksA) <- colnames(EntorhinalCortex_AD_H3K27ac_peaks)
  peaks <- data.frame(peaksA)
  counts <- data.frame(peaksB)
  AD_consensusPeaks <- ConsensusPeaks(peaks, counts,
                                      
                                      ref_H3K27ac_peaks_hg38, ref_H3K27ac_counts_hg38, bedtoolspath)
  
  AD_MF_noBAM <- CelltypeProportion(AD_consensusPeaks$newBulkTPM, AD_consensusPeaks$newRefTPM, AD_consensusPeaks$consensusPeaks,
                                    
                                    refSamples, NULL)
  setwd("~/Library/CloudStorage/OneDrive-UniversityofExeter/paper/CHAS/rerun/results/AlexiColours")
  pdf(paste0(cells[cell],"_CHAS_newColours_PW.pdf"), width = 20, height=15)
  
  
  # prepare dataset
  
  x <- ncol(AD_MF_noBAM[["proportions"]])-5
  
  names(AD_MF_noBAM[["proportions"]])[5+x] <- "Other"
  
  proportions <- AD_MF_noBAM[["proportions"]]*100
  
  # Plotting
  
  #par(mar=c(14.5, 5, 4, 9), xpd=TRUE)
 # par(mar=c(12, 5, 4, 8), xpd=TRUE)
   par(mar=c(16, 5, 4, 15), xpd=TRUE)
  
  par(las = 2)
  
  barplot(t(proportions),
          
          xlab = NULL, ylab = NULL, border = NA,
          cex.names = 1.6,
          
          col = c("#E69F00", "#BC3934" , "#009E73", "#0072B2",
                  
                  brewer.pal(8, "Accent")[0:x],"#DDDDDD"))
  
  legend("right",inset=c(-0.2,0),cex = 0.8,
         
         names(proportions),
         
         fill = c("#E69F00", "#BC3934" , "#009E73", "#0072B2",
                  
                  brewer.pal(8, "Accent")[0:x],"#DDDDDD"))
  
  #mtext("Samples", side=1, line=8, font=1, adj = 0.5, cex=1.2, las=1)
  
  mtext("Predicted proportions (%)", side=2, line=2.5,adj=0.5, font=1,cex=2.5, las=3)
  
  dev.off() 
  
  
  write.csv(proportions, file= paste0(cells[cell],"_celltypeproportions_PW.csv"), row.names = TRUE)
  
  ##  Cell type only plot ##
  AlexiColours <- c("OLIG2"="#0072B2","Astrocyte"= "#E69F00", "PU1"= "#BC3934", "NEUN"="#009E73", "otherCells"="#DDDDDD")

  

  
 # pdf(paste0(cells[cell],"_CHAS_CellOnly_PW.pdf"), width = 20, height=19)
  pdf(paste0(cells[cell],"_CHAS_CellOnly_PW.pdf"), width = 20, height=18)
  
  
  
  # prepare dataset
  
  x <- ncol(AD_MF_noBAM[["proportions"]])-5
  
  names(AD_MF_noBAM[["proportions"]])[5+x] <- "Other"
  
  proportions <- AD_MF_noBAM[["proportions"]]*100
  colnames(proportions) <- c("Astrocyte", "PU1", "NEUN", "OLIG2", "Other")  
  proportionsCELL <- proportions[colnames(proportions) %in%  cells[cell]]
  
  cellColour <- AlexiColours[which(names(AlexiColours)==cells[cell])]
  fullnames <- c("Microglia", "Oligodendrocytes", "Neurons")
  names(fullnames) <- cells
  name <- fullnames[which(names(fullnames)==cells[cell])]
  names(proportionsCELL)[1] <- name
  
  # Plotting
  
  par(mar=c(22, 5, 4, 8), xpd=TRUE)
  
  par(las = 2)
  
  
  barplot(t(proportionsCELL),
          ylim = c(0,100),
          xlab = NULL, ylab = NULL, border = NA,
          cex.names = 2.2,
          
          col = c(cellColour))
  
  legend("right",inset=c(-0.12,0),cex = 1,
         
         names(proportionsCELL),
         
         
         fill = c(cellColour))
  
  #mtext("Samples", side=1, line=8, font=1, adj = 0.5, cex=1.2, las=1)
  
  
  mtext("Predicted proportion (%)", side=2, line=2.5,adj=0.5, font=1,cex=2.5, las=3)
  
  dev.off() 

  ## CASE CONTROL COMPARISONS ##
  dfmCELL <- dfmm
  dfmCELL$cellID <- paste0(cells[cell],"_",dfmm$Donor,"_H3K27ac")
  dfmCELL$cellID <- gsub("NEUN", "NeuN", dfmCELL$cellID)
  dfmCELL$cellID <- gsub("OLIG2", "Olig2", dfmCELL$cellID)
  dfmAD <- dfmCELL[which(dfmCELL$AD==1),]
  dfmCONT <- dfmCELL[-which(dfmCELL$AD==1),]
  proportionsCELL$Cell <- cells[cell]
  proportionsCELL_AD <- proportionsCELL[which(row.names(proportionsCELL)%in%dfmAD$cellID),]
  proportionsCELL_AD$AD <- "AD"
  proportionsCELL_CONT<- proportionsCELL[which(row.names(proportionsCELL)%in%dfmCONT$cellID),]
  proportionsCELL_CONT$AD <- "CONT"
  celldf <-  paste0("dfviolin1_cells",cells[cell])
  celldf <- rbind(proportionsCELL_AD, proportionsCELL_CONT)
  colnames(celldf) <- c("Proportion", "Cell", "AD")
  dfviolin <- celldf
  write.csv(dfviolin, paste0("dfviolin_",cells[cell],".csv"))

  
}

## Violin, edited plot code from Yukyee ##
  data1 <- read.csv("dfviolin_PU1.csv")
  data2 <- read.csv("dfviolin_OLIG2.csv")
  data3 <- read.csv("dfviolin_NEUN.csv")
  data <- rbind(data1,data2,data3)
  head(data)
  data$Cell <- gsub("PU1", "Microglia", data$Cell)
  data$Cell <- gsub("OLIG2", "Oligodendrocytes", data$Cell)
  data$Cell <- gsub("NEUN", "Neurons", data$Cell)
  names(data)[4] <- "Group"
  data$Group <- gsub("CONT", "Control", data$Group)

  
violinPlot <- ggplot(data, aes(x = Group, y = Proportion, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), fill = "white") +
  geom_jitter(position = position_jitter(width = 0.15), size = 1, alpha = 0.7) +
  labs(title = "CHAS predicted cell type proportions within cell type enriched nuclei populations per group",
       x = "Disease Status",
       y = "Proportion (%) predicted cell type enrichment") +
  theme_minimal() +
  scale_fill_manual(values = c("AD" = "darkorange", "Control" = "darkgrey")) +
  facet_wrap(~ Cell, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5))

violinPlot
ggsave(plot = violinPlot, filename = "Violin_CHAS_groups_updatedPW.pdf", width = 12, height = 7)


# Wilcoxon test 
resultsWilcoxon_PU1 <- data1 %>%
  summarize(
    p_value = wilcox.test(Proportion ~ AD)$p.value
  )
#0.0501
resultsWilcoxon_OLIG2 <- data2 %>%
  summarize(
    p_value = wilcox.test(Proportion ~ AD)$p.value
  )
#0.37
resultsWilcoxon_NEUN <- data3 %>%
summarize(
  p_value = wilcox.test(Proportion ~ AD)$p.value
)
resultsWilcoxon <- rbind(resultsWilcoxon_NEUN,resultsWilcoxon_OLIG2, resultsWilcoxon_PU1)
row.names(resultsWilcoxon) <- c("NeuN", "Olig2", "PU1")
write.csv(resultsWilcoxon, "resultsWilcoxon.csv")





### NEUN case control##
write.csv(proportions_AD, "Proportions_NeuN_AD.csv")
write.csv(proportions_cont, "Proportions_NeuN_cont.csv")
## case control groups PU1 ##

## PU1 ##

dfmm <- read_csv("Metadata_per_file.csv")
colnames(dfc)
dfmm$AD <- ifelse(dfmm$`Braak   Stage`<3,0,1)
dfmAD <- dfmm[which(dfmm$AD==1),]
dfmCONT <- dfmm[-which(dfmm$AD==1),]
proportions <- read.csv("PU1_celltypeproportions_PW.csv")
row.names(proportions) <- proportions$X
proportions <- proportions[-1]
proportions_AD <- proportions[which(row.names(proportions)%in%dfmAD$NewID),]
proportions_CONT <- proportions[-which(row.names(proportions)%in%dfmAD$NewID),]
write.csv(proportions_AD, "proportions_PU1_AD.csv")
write.csv(proportions_CONT, "proportions_PU1_CONT.csv")
#  means
#54.11603  57.45774
## (higher in controls)


## Olig2 ##
proportions <- read.csv("OLIG2_celltypeproportions_PW.csv")
row.names(proportions) <- proportions$X
proportions <- proportions[-1]
proportions_AD <- proportions[which(row.names(proportions)%in%dfmAD$NewID),]
proportions_CONT <- proportions[-which(row.names(proportions)%in%dfmAD$NewID),]
write.csv(proportions_AD, "proportions_Olig2_AD.csv")
write.csv(proportions_CONT, "proportions_Olig2_CONT.csv")












