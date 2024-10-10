## plot volcanos from EdgeR results ###
setwd("~/Library/CloudStorage/OneDrive-UniversityofExeter/paper/EdgeR")
colours <- c("#E69F00", "#BC3934" , "#009E73", "#0072B2")
EdgeR_PU1 <- read_csv("~/Library/CloudStorage/OneDrive-UniversityofExeter/paper/EdgeR/EdgeR_PU1.csv")
EdgeR_NeuN <-read_csv("~/Library/CloudStorage/OneDrive-UniversityofExeter/paper/EdgeR/EdgeR_NeuN.csv")
EdgeR_Olig2 <-read_csv("~/Library/CloudStorage/OneDrive-UniversityofExeter/paper/EdgeR/EdgeR_Olig2.csv")


library(ggrepel)
library(ggplot2)
library(GenomicRanges)
# plot adding up all layers we have seen so far
# The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)

de <- EdgeR_PU1
de$Differential <- "Non-significant"
de$Differential[de$logFC < 0 & de$FDR < 0.05] <- "Differentially acetylated"
de$Differential[de$logFC > 0 & de$FDR < 0.05] <- "Differentially acetylated"
table(de$Differential)
colnames(de)
de <- as.data.frame(de)
de$Differential <-as.factor(de$Differential)

colnames(de)
library(ggplot2)
library(ggtext)
library(extrafont)
font_import()  # Please run this line to import the fonts; only need to do it once
loadfonts(quiet = T) 

pdf("Volcano_PU1_PW.pdf",, width=7.34, height =5.86)
### updated volcano ###
p <- ggplot(data=de, aes(x=logFC, y=-log10(FDR), col=Differential)) +
  geom_point(size=0.5) + 
  theme_minimal() +
  scale_color_manual(values=c("#BC3934" ,"mistyrose4", "#BC3934")) + 
  geom_hline(yintercept=-log10(0.05), col="black", linetype=2) +
  labs(color="") +
  theme( 
        legend.text = element_text(size=12))+
  labs(title = "<span style = 'font-size: 12pt'>**Differential acetylation:** Oligodendrocytes") + 
  theme(plot.title = element_textbox_simple(size = 10,
                                            #family = "Arial",
                                            color = "#007575",
                                            fill = "#e0f2f2",
                                            halign = 0.05,
                                            lineheight = 1.5,
                                            padding = margin(5, 1, 5, 1), 
                                            margin = margin(0, 0, 5, 0)))+
  theme(legend.position = "none")
 # guides(color = guide_legend(override.aes = list(size = 3) ) )
# labs(title = "<span style = 'font-size: 14pt'>Differential **A**cetylation **M**icroglia<i/></span><br>AD versus control PFC samples") + 
p
dev.off()
#####

## Olig2 ##
de <- EdgeR_Olig2
de$Differential <- "Non-significant"
de$Differential[de$logFC < 0 & de$FDR < 0.05] <- "Differentially acetylated"
de$Differential[de$logFC > 0 & de$FDR < 0.05] <- "Differentially acetylated"
table(de$Differential)
colnames(de)
de <- as.data.frame(de)
de$Differential <-as.factor(de$Differential)

colnames(de)
library(ggplot2)
library(ggtext)
library(extrafont)

pdf("Volcano_Olig2_PW.pdf",, width=7.34, height =5.86)
### updated volcano ###
p <- ggplot(data=de, aes(x=logFC, y=-log10(FDR), col=Differential)) +
  geom_point(size=0.5) + 
  theme_minimal() +
  scale_color_manual(values=c("#0072B2" ,"mistyrose4", "#0072B2")) + 
  geom_hline(yintercept=-log10(0.05), col="black", linetype=2) +
  labs(color="") +
  theme( 
    legend.text = element_text(size=12))+
  labs(title = "<span style = 'font-size: 12pt'>**Differential acetylation:** Oligodendrocytes") + 
  theme(plot.title = element_textbox_simple(size = 10,
                                            #family = "Arial",
                                            color = "#007575",
                                            fill = "#e0f2f2",
                                            halign = 0.05,
                                            lineheight = 1.5,
                                            padding = margin(5, 1, 5, 1), 
                                            margin = margin(0, 0, 5, 0)))+
  theme(legend.position = "none")
# guides(color = guide_legend(override.aes = list(size = 3) ) )
# labs(title = "<span style = 'font-size: 14pt'>Differential **A**cetylation **M**icroglia<i/></span><br>AD versus control PFC samples") + 
p
dev.off()

## NeuN #
## Olig2 ##
de <- EdgeR_NeuN
de$Differential <- "Non-significant"
de$Differential[de$logFC < 0 & de$FDR < 0.05] <- "Differentially acetylated"
de$Differential[de$logFC > 0 & de$FDR < 0.05] <- "Differentially acetylated"
table(de$Differential)
colnames(de)
de <- as.data.frame(de)
de$Differential <-as.factor(de$Differential)

colnames(de)
library(ggplot2)
library(ggtext)
library(extrafont)

pdf("Volcano_NeuN_PW.pdf",, width=7.34, height =5.86)
### updated volcano ###
p <- ggplot(data=de, aes(x=logFC, y=-log10(FDR), col=Differential)) +
  geom_point(size=0.5) + 
  theme_minimal() +
  scale_color_manual(values=c("#009E73" ,"mistyrose4", "#009E73")) + 
  geom_hline(yintercept=-log10(0.05), col="black", linetype=2) +
  labs(color="") +
  theme( 
    legend.text = element_text(size=12))+
  labs(title = "<span style = 'font-size: 12pt'>**Differential acetylation:** Neurons") + 
  theme(plot.title = element_textbox_simple(size = 10,
                                            #family = "Arial",
                                            color = "#007575",
                                            fill = "#e0f2f2",
                                            halign = 0.05,
                                            lineheight = 1.5,
                                            padding = margin(5, 1, 5, 1), 
                                            margin = margin(0, 0, 5, 0)))+
  theme(legend.position = "none")
# guides(color = guide_legend(override.aes = list(size = 3) ) )
# labs(title = "<span style = 'font-size: 14pt'>Differential **A**cetylation **M**icroglia<i/></span><br>AD versus control PFC samples") + 
p
dev.off()
