## ------------------------------------------------------------------------
library(ggplot2) 

## ------------------------------------------------------------------------
# Load the data 
data(iris)

# Size of the dataframe 
dim(iris)

# First few lines of the dataset 
head(iris) 

## ------------------------------------------------------------------------
# Summary statistics 
summary(iris) 

## ------------------------------------------------------------------------
pca_iris <- prcomp(iris[,-5]) 

## ------------------------------------------------------------------------
summary(pca_iris) 

## ------------------------------------------------------------------------
plot_iris <- data.frame(pca_iris$x, Species=iris$Species) 
head(plot_iris)

## ----fig.width=7.5, fig.height=4.5---------------------------------------
ggplot(plot_iris, aes(x=PC1, y=PC2, color=Species)) + geom_point() 

## ----fig.width=7.5, fig.height=4.5---------------------------------------
ggplot(plot_iris, aes(x=PC3, y=PC4, color=Species)) + geom_point() 

## ----fig.width=7.5, fig.height=4.5---------------------------------------
with(plot_iris, plot(x=PC1, y=PC2, col=Species)); grid() 

## ------------------------------------------------------------------------
# Load the data 
library(pdfCluster) 
data(oliveoil)

# Size of the dataframe 
dim(oliveoil)

# First few lines of the dataset 
head(oliveoil) 

## ------------------------------------------------------------------------
# Summary statistics 
summary(oliveoil) 

## ------------------------------------------------------------------------
pca_oil <- prcomp(oliveoil[,-c(1,2)]) 
summary(pca_oil) 

## ------------------------------------------------------------------------
plot_oil <- data.frame(pca_oil$x, oliveoil[,1:2]) 
head(plot_oil) 

## ----fig.width=7.5, fig.height=4.5---------------------------------------
ggplot(plot_oil, aes(x=PC1, y=PC2, color=region, shape=macro.area)) + geom_point() + coord_fixed() 

## ----fig.width=7.5, fig.height=4.5---------------------------------------
ggplot(plot_oil, aes(x=PC3, y=PC4, color=region, shape=macro.area)) + geom_point() + coord_fixed() 

## ------------------------------------------------------------------------
library(Shenzen2016)
data("zebrafish")

## ------------------------------------------------------------------------
above_threshold <- rowSums(zebrafish$Expression >= 5) >= 3
EM_zebra <- subset(zebrafish$Expression, above_threshold)

## ------------------------------------------------------------------------
library(edgeR)

# Create dge list object
dge_zebra <- DGEList(EM_zebra)

# Calculate normalization factors
dge_zebra <- calcNormFactors(dge_zebra, method="TMM")

# Normalize the expression values
logTMM_zebra <- cpm(dge_zebra, log=TRUE, prior.count=3)

## ------------------------------------------------------------------------
pca_zebra <- prcomp(t(logTMM_zebra))
summary(pca_zebra)

## ------------------------------------------------------------------------
plot_zebra <- data.frame(pca_zebra$x, zebrafish$Design) 
head(plot_oil) 

## ----fig.width=7.5, fig.height=4.5---------------------------------------
ggplot(plot_zebra, aes(x=PC1, y=PC2, color=gallein)) + geom_point() + coord_fixed() 

## ----fig.width=7.5, fig.height=4.5---------------------------------------
ggplot(plot_zebra, aes(x=PC3, y=PC4, color=gallein)) + geom_point() + coord_fixed()

## ------------------------------------------------------------------------
# Load the data
data("yeast")

# Trim 
EM_yeast <- subset(yeast$Expression, rowSums(yeast$Expression >= 5) >= 3)

# Calculate normalization factors
dge_yeast <- calcNormFactors(DGEList(EM_yeast), method="TMM")

# Normalize the expression values
logTMM_yeast <- cpm(dge_yeast, log=TRUE, prior.count=3)

# PCA
pca_yeast <- prcomp(t(logTMM_yeast))
plot_yeast <- data.frame(pca_yeast$x, yeast$Design) 

## ----fig.width=7.5, fig.height=4.5---------------------------------------
ggplot(plot_yeast, aes(x=PC1, y=PC3, color=minute, shape=strain)) + geom_point() + coord_fixed() 

## ------------------------------------------------------------------------
# Load the data
data("pasilla")

# Trim 
EM_pasilla <- subset(pasilla$Expression, rowSums(pasilla$Expression >= 5) >= 3)

# Calculate normalization factors
dge_pasilla <- calcNormFactors(DGEList(EM_pasilla), method="TMM")

# Normalize the expression values
logTMM_pasilla <- cpm(dge_pasilla, log=TRUE, prior.count=3)

# PCA
pca_pasilla <- prcomp(t(logTMM_pasilla))
plot_pasilla <- data.frame(pca_pasilla$x, pasilla$Design) 

## ----fig.width=7.5, fig.height=4.5---------------------------------------
ggplot(plot_pasilla, aes(x=PC1, y=PC2, color=condition, shape=type)) + geom_point() + coord_fixed() 

