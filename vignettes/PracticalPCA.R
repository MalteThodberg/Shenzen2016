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
ggplot(plot_iris, aes(x=PC1, y=PC2, color=Species)) + geom_point()  + coord_fixed()

## ----fig.width=7.5, fig.height=4.5---------------------------------------
ggplot(plot_iris, aes(x=PC3, y=PC4, color=Species)) + geom_point() + coord_fixed()

## ----fig.width=7.5, fig.height=4.5---------------------------------------
with(plot_iris, plot(x=PC1, y=PC2, col=Species)); grid() 

## ------------------------------------------------------------------------
# Load the data 
suppressPackageStartupMessages(library(pdfCluster))
data(oliveoil)

# Size of the dataframe 
dim(oliveoil)

# First few lines of the dataset 
head(oliveoil) 

## ------------------------------------------------------------------------
# Summary statistics 
summary(oliveoil) 

## ------------------------------------------------------------------------
# PCA without the first two columns
pca_oil <- prcomp(oliveoil[,-c(1,2)]) 

# Save the data for plotting
plot_oil <- data.frame(pca_oil$x, oliveoil[,1:2]) 

## ----fig.width=7.5, fig.height=4.5---------------------------------------
ggplot(plot_oil, aes(x=PC1, y=PC2, color=region, shape=macro.area)) + geom_point() + coord_fixed() 

## ------------------------------------------------------------------------
library(Shenzen2016)
data("zebrafish")

## ------------------------------------------------------------------------
above_threshold <- rowSums(zebrafish$Expression >= 10) >= 3
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
data("pasilla")
data("tissues")

# Inspect the dimensionality
dim(yeast$Expression)

# Inspect the design by
yeast$Design


## ------------------------------------------------------------------------
# Set a dataset
dataset <- tissues

# Trim: Play around with numbers here 
EM <- subset(dataset$Expression, rowSums(dataset$Expression >= 5) >= 3)

# Calculate normalization factors: Play around with the method argument
dge <- calcNormFactors(DGEList(EM), method="TMM")

# Normalize the expression values
logTMM <- cpm(dge, log=TRUE, prior.count=5)

# Perform the PCA and save data for plotting
pca <- prcomp(t(logTMM))
plot_data <- data.frame(pca$x, dataset$Design) 

## ------------------------------------------------------------------------
ggplot(plot_data, aes(x=PC1, y=PC2, color=cell.type, label=cell.type)) + geom_point() + coord_fixed() 

