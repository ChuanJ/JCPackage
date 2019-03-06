# JCPackage
Adjust for positional and batch effects using ComBat function

In this package, we will carefully compare the effects of positions and batches on your datasets. Then a appropriate correction order will be provided.

You should provide the sample Sentrix information and the methylation levels in matrix file without NA file. The example files shown in the Data file.

## Install from GitHub:
```R
library(devtools)
install_github("ChuanJ/JCpackage")

library(posibatch)
```
## Usage:
### input the methylation levels
```R
dat <-as.matrix(read.csv("beta.csv",header=T,sep=",",row.names=1))
```
### input the Sentrix information
```R
Sentrix <- read.csv("Sentrix.csv",header=T,sep=",")
```
### Correct the positional effect and batch effect
```R
results <- posibatches(dat, Sentrix, posi=TRUE, batch=TRUE, par.prior=TRUE, prior.plots=FALSE, mean.only.posi=FALSE, mean.only.batch=FALSE)
```
### output the file
```R
write.csv(results,"dataAfterBatch&position.csv")
```
