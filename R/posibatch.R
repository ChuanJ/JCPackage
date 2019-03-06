#' Adjust for positional effects and batch effects using posibatch in Illumina Beadchips
#' 
#' posibatch allows users to adjust for positional effects and batch effects in Illumina Beadchips. Positional effects means the same sample in
#'  different physical positions on the array could be measured as different methylation or expression levels, and batch effects. Batch effects
#'  are sub-groups of measurements that have qualitatively different behaviour across conditions and are unrelated to the biological or
#'  scientific variables in a study. It uses either parametric or non-parametric empirical Bayes frameworks for adjusting data for positional
#'  effects and batch effects. Users are returned an expression or methylation matrix that has been corrected for position effects and batch effects. 
#'  The input data are assumed to be cleaned and normalized before batch effect removal. 
#' 
#' @param dat Genomic measure matrix (dimensions probe x sample) - for example, expression matrix
#' @param Sentrix The list of position numbers and chip numbers - for example, sampleID, sampleNames (8963303124_R01C01) and batches
#' @param posi {Position covariate}
#' @param batch {Batch covariate}
#' @param mod Model matrix for outcome of interest and other covariates besides batch and position
#' @param par.prior (Optional) TRUE indicates parametric adjustments will be used, FALSE indicates non-parametric adjustments will be used
#' @param prior.plots (Optional) TRUE give prior plots with black as a kernel estimate of the empirical batch effect density and red as the parametric
#' @param mean.only.posi (Optional) FALSE If TRUE posibatch only corrects the mean of the positional effect (no scale adjustment). If one position has only one sample, setting mean.only=TRUE.
#' @param mean.only.batch (Optional) FALSE If TRUE posibatch only corrects the mean of the batch effect (no scale adjustment). If one batch has only one sample, setting mean.only=TRUE.
#'
#' @return data A probe x sample genomic measure matrix, adjusted for positional effects and batch effects.
#' 
#' @export
#' 


################################################################################
################################################################################
posibatches <- function(dat, Sentrix, posi=TRUE, batch=TRUE, par.prior=TRUE, prior.plots=FALSE, mean.only.posi=FALSE, mean.only.batch=FALSE) {
  require(sva)
  require(lme4)
  #get the position and batch information----------------------------------------
  ### Extraction of the position numbers and chip numbers--------------------
  #if (is.null(SentrixVector)){
  if (is.null(Sentrix)){
    stop('Sentrix informations must be provided.')
  }
  #chips<-as.numeric(factor(substr(Sentrix$sampleNames, 1, 10))) 
  #positions<-as.numeric(factor(substr(Sentrix$sampleNames, 12, 17))) 
  
  chips <- as.numeric(factor(sapply(strsplit(as.character(Sentrix$sampleNames), "_"), function(x) x[1])))
  positions <- as.numeric(factor(sapply(strsplit(as.character(Sentrix$sampleNames), "_"), function(x) x[2])))
                               
  if (length(positions)!=length(chips)){
    stop('positions and chips must have the same length')
  }
  if (sum(positions>12)>0){
    stop('Position number cannot be greater than 12')
  }
  if (sum(is.na(chips))>1 || sum(is.na(positions))>1){ 
    stop('One or more position or chip numbers missing')
  }
  
  if (sum(is.na(Sentrix$batches))==length(chips)){
    batches<-chips
  } else {
    batches<-as.factor(Sentrix$batches)
  }
  ################################################################################
  ################################################################################
  pct_threshold = .8 # Amount of variability desired to be explained by the principal components.  Set to match the results in book chapter and SAS code.  User can adjust this to a higher (>= 0.8) number but < 1.0
  dataRowN <- nrow(dat)
  dataColN <- ncol(dat)
  
  ########## Center the data (center rows) ##########
  datCentered <- matrix(data = 0, nrow = dataRowN, ncol = dataColN)
  datCentered_transposed = apply(dat, 1, scale, center = TRUE, scale = FALSE)
  datCentered = t(datCentered_transposed)
  
  exp_design<-data.frame(cbind(positions,batches))
  expDesignRowN <- nrow(exp_design)
  expDesignColN <- ncol(exp_design)
  myColNames <- names(exp_design)
  
  
  ########## Compute correlation matrix ##########
  
  theDataCor <- cor(datCentered)
  
  ########## Obtain eigenvalues ##########
  
  eigenData <- eigen(theDataCor)
  eigenValues = eigenData$values
  ev_n <- length(eigenValues)
  eigenVectorsMatrix = eigenData$vectors
  eigenValuesSum = sum(eigenValues)
  percents_PCs = eigenValues /eigenValuesSum 
  
  ########## Merge experimental file and eigenvectors for n components ##########
  
  my_counter_2 = 0
  my_sum_2 = 1
  for (i in ev_n:1){
    my_sum_2  = my_sum_2 - percents_PCs[i]
    if ((my_sum_2) <= pct_threshold ){
      my_counter_2 = my_counter_2 + 1
    }
    
  }
  if (my_counter_2 < 3){
    pc_n  = 3
    
  }else {
    pc_n = my_counter_2 
  }
  
  # pc_n is the number of principal components to model
  
  pc_data_matrix <- matrix(data = 0, nrow = (expDesignRowN*pc_n), ncol = 1)
  mycounter = 0
  for (i in 1:pc_n){
    for (j in 1:expDesignRowN){
      mycounter <- mycounter + 1
      pc_data_matrix[mycounter,1] = eigenVectorsMatrix[j,i]
      
    }
  }
  
  AAA <- exp_design[rep(1:expDesignRowN,pc_n),]
  
  Data <- cbind(AAA,pc_data_matrix)
  
  ####### Edit these variables according to your factors #######
  
  variables <- c(colnames(exp_design))
  for (i in 1:length(variables)) {
    Data$variables[i] <- as.factor(Data$variables[i])
  }
  
  ########## Mixed linear model ##########
  op <- options(warn = (-1))
  #effects_n = expDesignColN + choose(expDesignColN, 2) + 1
  effects_n = expDesignColN  + 1
  randomEffectsMatrix <- matrix(data = 0, nrow = pc_n, ncol = effects_n)
  
  model.func <- c()
  index <- 1
  for (i in 1:length(variables)) {
    mod = paste("(1|", variables[i], ")", sep = "")
    model.func[index] = mod
    index = index + 1
  }
  
  function.mods <- paste(model.func, collapse = " + ")
  
  for (i in 1:pc_n) {
    y = (((i - 1) * expDesignRowN) + 1)
    funct <- paste("pc_data_matrix", function.mods, sep = " ~ ")
    Rm1ML <- lmer(funct, Data[y:(((i - 1) * expDesignRowN) + 
                                   expDesignRowN), ], REML = TRUE, control=lmerControl(check.nobs.vs.nlev = "ignore",check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"),verbose = FALSE, 
                  na.action = na.omit)
    randomEffects <- Rm1ML
    randomEffectsMatrix[i, ] <- c(unlist(VarCorr(Rm1ML)), 
                                  resid = sigma(Rm1ML)^2)
  }
  effectsNames <- c(names(getME(Rm1ML, "cnms")), "resid")
  ########## Standardize Variance ##########
  
  randomEffectsMatrixStdze <- matrix(data = 0, nrow = pc_n, ncol = effects_n)
  for (i in 1:pc_n){
    mySum = sum(randomEffectsMatrix[i,])
    for (j in 1:effects_n){
      randomEffectsMatrixStdze[i,j] = randomEffectsMatrix[i,j]/mySum	
    }
  }
  
  ########## Compute Weighted Proportions ##########
  
  randomEffectsMatrixWtProp <- matrix(data = 0, nrow = pc_n, ncol = effects_n)
  for (i in 1:pc_n){
    weight = eigenValues[i]/eigenValuesSum
    for (j in 1:effects_n){
      randomEffectsMatrixWtProp[i,j] = randomEffectsMatrixStdze[i,j]*weight
    }
  }
  ######### Compute Weighted Ave Proportions ##########
  
  randomEffectsSums <- matrix(data = 0, nrow = 1, ncol = effects_n)
  randomEffectsSums <-colSums(randomEffectsMatrixWtProp)
  totalSum = sum(randomEffectsSums)
  randomEffectsMatrixWtAveProp <- matrix(data = 0, nrow = 1, ncol = effects_n)
  
  for (j in 1:effects_n){
    randomEffectsMatrixWtAveProp[j] = randomEffectsSums[j]/totalSum 	
    
  }
  
  if(randomEffectsMatrixWtAveProp[,1]<randomEffectsMatrixWtAveProp[,2]){
    if(posi==TRUE){
      afterposiExp<-ComBat(dat = dat, batch = positions, mod = mod, par.prior = par.prior, prior.plots = prior.plots, mean.only = mean.only.posi)
      dat<-afterposiExp
    }
    if(batch==TRUE){
      afterbatchExp<-ComBat(dat = dat, batch = batches, mod = mod, par.prior = par.prior, prior.plots = prior.plots, mean.only = mean.only.batch)
      dat<-afterbatchExp
    }
    else{
      if(batch==TRUE){
        afterbatchExp<-ComBat(dat = dat, batch = batches, mod = mod, par.prior = par.prior, prior.plots = prior.plots, mean.only = mean.only.batch)
        dat<-afterbatchExp
      }
      if(posi==TRUE){
        afterposiExp<-ComBat(dat = dat, batch = positions, mod = mod, par.prior = par.prior, prior.plots = prior.plots, mean.only = mean.only.posi)
        dat<-afterposiExp
      }
    }
     
return(dat)
  }
  
  ################################################################################
  ################################################################################

  
}







