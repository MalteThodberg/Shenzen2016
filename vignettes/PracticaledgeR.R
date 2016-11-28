## ------------------------------------------------------------------------
dummy <- data.frame(factor1=c("A", "A", "B", "B", "A", "A", "B", "B"),
										 factor2=c("C", "C", "C", "C", "D", "D", "D", "D"))

dummy

## ------------------------------------------------------------------------
model.matrix(~factor1, data=dummy)

## ------------------------------------------------------------------------
model.matrix(~factor1+factor2, data=dummy)

## ------------------------------------------------------------------------
model.matrix(~factor1+factor2+factor1:factor2, data=dummy)

## ------------------------------------------------------------------------
model.matrix(~factor1*factor2, data=dummy)

## ------------------------------------------------------------------------
# Load the data
library(Shenzen2016)
library(ggplot2)
library(edgeR)
data("zebrafish")

# Trim 
EM_zebra <- subset(zebrafish$Expression, rowSums(zebrafish$Expression >= 5) >= 3)

# Calculate normalization factors
dge_zebra <- DGEList(EM_zebra)
dge_zebra <- calcNormFactors(dge_zebra, method="TMM")

## ------------------------------------------------------------------------
zebrafish$Design

## ------------------------------------------------------------------------
mod <- model.matrix(~gallein, data=zebrafish$Design)

## ------------------------------------------------------------------------
disp_zebra <- estimateDisp(dge_zebra, design=mod, robust=TRUE)

## ------------------------------------------------------------------------
plotBCV(y=disp_zebra)

## ------------------------------------------------------------------------
fit_zebra <- glmFit(disp_zebra, design=mod)

## ------------------------------------------------------------------------
# Using the name of the coefficient
gallein <- glmLRT(fit_zebra, coef="galleintreated")

# Using the index of the coefficient
gallein <- glmLRT(fit_zebra, coef=2)

## ------------------------------------------------------------------------
topTags(gallein)

## ------------------------------------------------------------------------
is_de <- decideTestsDGE(gallein, p.value=0.1)
head(is_de)
summary(is_de)

## ------------------------------------------------------------------------
plotSmear(gallein, de.tags=rownames(gallein)[is_de != 0])

## ------------------------------------------------------------------------
res <- as.data.frame(topTags(gallein, n=Inf, sort.by="none"))
ggplot(res, aes(x=logFC, y=-log10(PValue), color=FDR < 0.05)) + geom_point()

## ------------------------------------------------------------------------
fitQL <- glmQLFit(y=disp_zebra, design=mod, robust=TRUE)
galleinQL <- glmQLFTest(fitQL, coef="galleintreated")
summary(decideTestsDGE(galleinQL, p.value=0.1))

## ------------------------------------------------------------------------
plotQLDisp(fitQL)

## ------------------------------------------------------------------------
# Load the data
data("yeast")
data("pasilla")

# Set a datasets
dataset <- pasilla

# Experiment with model matrices
mod <- model.matrix(~condition, data=dataset$Design)

# Look how the experiment is designed
dataset$Design

## ------------------------------------------------------------------------
# Trim 
EM <- subset(dataset$Expression, rowSums(dataset$Expression >= 10) >= 3)

# Calculate normalization factors
dge <- calcNormFactors(DGEList(EM), method="TMM")

# Calculate dispersion
disp <- estimateDisp(dge, design=mod, robust=TRUE)

# Fit models
fit <- glmFit(disp, design=mod)

# Perform the test of the given coefficient
res <- glmLRT(fit, coef="conditionuntreated")

## ---- eval=FALSE, include=TRUE-------------------------------------------
#  # Dispersion
#  plotBCV(disp)
#  
#  # topTages
#  topTags(res)
#  
#  # Up and down regulated genes
#  summary(decideTestsDGE(res))
#  
#  # Smear plot
#  plotSmear(res, de.tags=rownames(res)[decideTestsDGE(res) != 0])
#  
#  # Volcano
#  all_genes <- as.data.frame(topTags(res, n=Inf, sort.by="none"))
#  ggplot(all_genes, aes(x=logFC, y=-log10(PValue), color=FDR < 0.05)) + geom_point()

