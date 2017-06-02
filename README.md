# JCPackage
Adjust for positional and batch effects using ComBat

You should provide the sample Sentrix information and the methylation levels in matrix file without NA file. The example files shown in the Data file.

## Install from GitHub:
```R
library(devtools)
install_github("ChuanJ/JCpackage")
```
## Usage:
### input the methylation levels
```R
dat <-as.matrix(read.csv("beta.csv",header=T,sep=",",row.names=1))
```
### input the Sentrix information
```R
Sentrix <- read.csv("Sentrix.csv",header=T,sep=",")
results <- posibatches(dat, Sentrix, posi=TRUE, batch=TRUE, par.prior=TRUE, prior.plots=FALSE, mean.only.posi=FALSE, mean.only.batch=FALSE)
```
### output the file
```R
write.csv(results,"dataAfterBatch&position.csv")
```
