## "dpIEL-Microbiome project, Batch correction"
### "Puspendu Sardar, Ph.D, Department of medicine, University of Cambridge, UK"

# install.packages("ggplot2")
library(ggplot2)

# Draw divergent plot
# Color based on value
color2 <- ifelse(divergent_data$Values < 0, "red", "blue")

pdf(file = "divergentPlot.pdf", width = 15.0, height = 7.0)
ggplot(divergent_data, aes(x = reorder(Species, Values), y = Values)) +
  geom_bar(stat = "identity",
           show.legend = FALSE,
           fill = color2,     
           color = "white") +
  geom_hline(yintercept = 0, color = 1, lwd = 0.2) +
  geom_text(aes(label = Species, # Text with groups
                hjust = ifelse(Values < 0, 1.0, 0),
                vjust = 0.5), size = 3.0) +
  xlab("Species") +
  ylab("Mean difference") +
  scale_y_continuous(breaks = seq(-6.5, 8.0, by = 0.8),
                     limits = c(-6.5, 8.0)) +
  coord_flip() +
  theme_minimal() +
  theme(axis.text.y = element_blank(),  # Remove Y-axis texts
        axis.ticks.y = element_blank(), # Remove Y-axis ticks
        panel.grid.major.y = element_blank(), # Remove horizontal grid
        panel.grid.minor.x = element_blank(), # Remove vertical grid
        panel.grid.major.x = element_blank()) # Remove vertical grid
dev.off()


# Draw Ellipses and Star plot
library(factoextra)

pca <- prcomp(Propionate_ordination, scale.=T, retx=T)
scores <-data.frame(pca$x)
write.table(scores, file = "Score_PCA_Propionate_ORF_level.tsv", quote = FALSE, sep = "\t", dec = ".")

# Ellipses plot
pdf(file = "Propionate_ECs.pdf", width = 10, height = 6)
fviz_pca_ind(pca, label="none", habillage=metadata_HBI$HBI_3Low,
             addEllipses=TRUE, ellipse.level=0.95, ellipse.type = "confidence",
             geom = "text", legend.title = "HBI") +
  #labs(title ="PCA Propionate biosynthesis", x = "PC1 (60.8%)", y = "PC2 (32.4%)") +
  theme_bw() +
  geom_point(aes(color=metadata_HBI$HBI_3Low), size=2)
dev.off()

# Star plot
#build ggplot dataframe with points (x,y) and corresponding groups (cluster)
gg <- data.frame(HBI=factor(metadata_HBI$HBI_3Low), x=scores$PC1, y=scores$PC2)

# calculate group centroid locations
centroids <- aggregate(cbind(x,y) ~ HBI, data=gg, mean)

# merge centroid locations into ggplot dataframe
gg <- merge(gg, centroids, by="HBI",suffixes=c("",".centroid"))
# generate star plot...

pdf(file = "PCA_star_plot_Propionate_pathway.pdf", width = 9, height = 6)
ggplot(gg) +
  geom_point(aes(x=x, y=y, color=HBI), size=2) +
  geom_point(data=centroids, aes(x=x, y=y, color=HBI), size=4) +
  geom_segment(aes(x=x.centroid, y=y.centroid, xend=x, yend=y, color=HBI)) +
  labs(title ="Propionate pathway", x = "PC1 (58.0%)", y = "PC2 (33.4%)") +    # The axes % values should be available
  theme_bw()
dev.off()

#### t-test on the axes
scores <- data.frame(res.pca$x)
scores$HBI_3Low <- metadata_HBI$HBI_3Low
##### Two-tail t-test
t.test(PC1 ~ HBI_3Low, data = scores)
##### One-tail t-test (alternative hypothesis)
t.test(PC1 ~ HBI_3Low, data = scores, alternative = "l")

boxplot(PC1 ~ HBI_3Low, data = scores, boxwex = 0.1, col = c("brown1", "cadetblue3"))