## "dpIEL-Microbiome project, Plotting human transcriptome data"
### "Puspendu Sardar, Ph.D, Department of medicine, University of Cambridge, UK"

plot1 <- ggplot(NULL, aes(x=sample_info$HBI_3Low,
                        y=log2(t(norm_counts["A2M",]+1)))) +
  geom_jitter(aes(shape=sample_info$HBI_3Low,
                  color=sample_info$HBI_3Low), size=3)+
  xlab(NULL) +
  ylab("Total expression \n log2 (norm counts +1)") +
  theme(legend.position = "bottom") +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title =element_text(size = 25),
        legend.position = 'none') +
  stat_summary(fun=mean,
               geom="point",
               shape= '_',
               size=14,
               colour= c('#b53d35', '#066e70'))


#### Plot PCA-biplot
library(factoextra)
pca_dat <- t(log2(norm_counts[c("CD28","CD44","CD58","CD69","CD83","CD86","CD8A","CTSW","IFNG","IL2RB","ITGAE","TNF"),]+1))
res.pca <- prcomp(pca_dat, scale = TRUE)

pdf(file = "deseq_norm_curated_genes_pca.pdf", width = 9, height = 6)
pca_plot <- fviz_pca_biplot(res.pca, col.ind = sample_info$HBI_3Low, geom = "point", addEllipses = TRUE, repel = TRUE,
                      ellipse.level=0.99, legend.title = "Disease score", ellipse.type = "confidence",
                     title='Principal Component Analysis of the curate genes') +
  theme_bw() +
  geom_point(aes(color=sample_info$HBI_3Low, shape=sample_info$HBI_3Low), size=1.5)
dev.off()

pca_plot

fviz_pca_ind(res.pca, label="none", habillage=propionate_ordi$HBI_3Low,
             addEllipses=TRUE, ellipse.level=0.9999, ellipse.type = "confidence",
             geom = "text", legend.title = "Disease score") +
  #labs(title ="PCA Core acetate pathway", x = "PC1 (63.6%)", y = "PC2 (22.6%)") +
  theme_bw() +
  geom_point(aes(color=propionate_ordi$HBI_3Low), size=2)

#### t-test on the axes
scores <- data.frame(res.pca$x)
scores$HBI_3Low <- metadata_HBI$HBI_3Low
##### Two-tail t-test
t.test(PC1 ~ HBI_3Low, data = scores)
##### One-tail t-test (alternative hypothesis)
t.test(PC1 ~ HBI_3Low, data = scores, alternative = "l")

boxplot(PC1 ~ HBI_3Low, data = scores, boxwex = 0.1, col = c("brown1", "cadetblue3"))


#### 3-D plot of the selected genes
library(canvasXpress)
data=t(log2(norm_counts[c("CD8A","ITGAE","CTSW"),]+1))

canvasXpress(data,
varAnnot=as.data.frame(sample_info$HBI_3Low, row.names = rownames(sample_info)),
axisTickScaleFontFactor=0.6,
axisTitleScaleFontFactor=0.6,
ellipseBy="sample_info$HBI_3Low",
colorBy="sample_info$HBI_3Low",
colorKey=list("sample_info$HBI_3Low"=list("High"="#F8766D", "Low"="#00BFC4")),
graphType="Scatter3D",
title="3D scatter plot of selected genes",
xAxis=list("CD8A"),
yAxis=list("ITGAE"),
zAxis=list("CTSW"),
showLoessFit = FALSE)