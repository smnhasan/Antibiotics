###################Water Quality in Sylhet City Corporation Area####################
#                            Mohammad Nayeem Hasan                                 #
####################################################################################


library(extrafont)
library(ggplot2)
library(pastecs)
library(corrplot)
library(ppcor)
library(factoextra)
library(psych)
library(GPArotation)
library(Hmisc)
library(dplyr)
library(ape)
library(plotly)

#loding main data
dat <- read.table("E:\\Bio-informatics\\PCA_CA - Antibiotic\\DataNanalysis\\data.csv",sep=',',header=T)
str(dat)

#descriptive statistics
stat.desc(dat)


#Box plot

#main data for boxplot
box <- read.table("E:\\Bio-informatics\\PCA_CA - Antibiotic\\DataNanalysis\\boxplot.csv",sep=',',header=T)
str(box)

p1 <- ggplot(box, aes(x = ï..Variables, y = Values)) +
  geom_boxplot(colour = "black", fill = "#56B4E9") +
  scale_y_continuous(name = "Values",
                     breaks = seq(0, 50, 10),
                     limits=c(0, 50)) +
  scale_x_discrete(name = "Antibiotic") +
  theme(axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title=element_text(size = 12, family="xkcd-Regular"),
        text=element_text(size = 12, family="xkcd-Regular"),
        axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 10))
p1

#Spearman Correlation Coefficient among Hydrological Parameters
antibiotics <- data.frame(dat$AMP, dat$PIT, dat$CTX, dat$CTR, dat$GEN,
                       dat$CX, dat$AZM, dat$TE, dat$COT, dat$CPM, dat$IMP,
                       dat$LE, dat$CXM)
colnames(antibiotics) <- c("AMP", "PIT","CTX", "CTR","GEN", "CX","AZM", "TE",
                        "COT", "CPM", "IMP", "LE", "CXM")
cor_sp = rcorr(as.matrix(antibiotics), type=c("spearman"))
cor_sp_coeff = cor_sp$r
round(cor_sp_coeff,2)
cor_sp_pval = cor_sp$P
round(cor_sp_pval,2)


#partial correlation All together
round(pcor(antibiotics, method = c("spearman"))$estimate,2)
round(pcor(antibiotics, method = c("spearman"))$p.value,2)

#Principal conponents analysis
res.pca <- prcomp(antibiotics, scale = T, center = T)

fviz_eig(res.pca,main = "Screeplot of the 10 PCs",
         addlabels = TRUE, 
         ylim = c(0, 30))

eigs <- res.pca$sdev^2
eigs
summary(res.pca)
pca.loadings <- res.pca$rotation
round(pca.loadings,3)

p <- plot_ly(
  type = 'scatterpolar',
  fill = "toself",
  opacity = 0.8
) %>%
  add_trace(
    r = c(0.312,0.234,0.384,0.417,0.283,0.310,0.288,0.079,-0.156,0.366,0.271,0.132,0.109),
    theta = c('AMP','PIT','CTX','CTR','GEN','CX','AZM','TE','COT','CPM','IMP','LE','CXM'),
    name = 'PC1 3.70'
  ) %>%
  add_trace(
    r = c(-0.119,0.126,0.198,0.118,-0.286,-0.101,-0.177,-0.574,-0.379,0.073,0.269,-0.492,-0.026),
    theta = c('AMP','PIT','CTX','CTR','GEN','CX','AZM','TE','COT','CPM','IMP','LE','CXM'),
    name = 'PC2 1.92'
  ) %>%
  add_trace(
    r = c(0.132,0.421,-0.077,-0.350,0.017,-0.024,0.266,0.283,-0.259,-0.328,0.323,-0.117,-0.482),
    theta = c('AMP','PIT','CTX','CTR','GEN','CX','AZM','TE','COT','CPM','IMP','LE','CXM'),
    name = 'PC3 1.44'
  ) %>%
  add_trace(
    r = c(-0.524,0.070,-0.090,0.114,0.438,-0.532,0.408,-0.030,0.055,0.202,0.041,-0.111,-0.013),
    theta = c('AMP','PIT','CTX','CTR','GEN','CX','AZM','TE','COT','CPM','IMP','LE','CXM'),
    name = 'PC4 1.36'
  ) %>%
  add_trace(
    r = c(0.249,-0.495,-0.016,0.013,0.037,-0.060,0.183,-0.122,-0.309,0.293,-0.327,0.011,-0.592),
    theta = c('AMP','PIT','CTX','CTR','GEN','CX','AZM','TE','COT','CPM','IMP','LE','CXM'),
    name = 'PC5 1.22'
  ) %>%
  layout(
    polar = list(
      radialaxis = list(
        visible = T,
        range = c(-0.5,0.5)
      )
    )
  )
p



#cluster Analysis

antibiotic_sta <- data.frame(dat$ï..ID, dat$AMP, dat$PIT, dat$CTX, dat$CTR, dat$GEN,
                           dat$CX, dat$AZM, dat$TE, dat$COT, dat$CPM, dat$IMP,
                           dat$LE, dat$CXM)
colnames(antibiotic_sta) <- c("ID","AMP", "PIT","CTX", "CTR","GEN", "CX","AZM", "TE",
                            "COT", "CPM", "IMP", "LE", "CXM")
antibiotic_sta <- t(scale(antibiotic_sta[2:14]))

str(antibiotic_sta)
head(antibiotic_sta)

# Compute dissimilarity matrix
res.dist <- dist(antibiotic_sta, method = "euclidean")

# Visualize the dissimilarity matrix
fviz_dist(res.dist, lab_size = 8)

# Enhanced hierarchical clustering
res.hc <- eclust(antibiotic_sta, "hclust" , k=5) # compute hclust
fviz_dend(res.hc, cex = 1, 
          main = "Dendogram Using Ward Linkage\nRescaled Distance",
          xlab = "Sampling Sites", ylab = "Distance", sub = "",rect = TRUE)

fviz_silhouette(res.hc)
fviz_cluster(res.hc) # scatter plot
