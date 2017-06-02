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
  #get the position and batch information----------------------------------------
  ### Extraction of the position numbers and chip numbers--------------------
  #if (is.null(SentrixVector)){
  if (is.null(Sentrix)){
    stop('Sentrix informations must be provided.')
  }
  chips<-as.numeric(factor(substr(Sentrix$sampleNames, 1, 10))) 
  positions<-as.numeric(factor(substr(Sentrix$sampleNames, 12, 17))) 
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
  if(posi==TRUE){
    afterposiExp<-ComBat(dat = dat, batch = positions, mod = mod, par.prior = par.prior, prior.plots = prior.plots, mean.only = mean.only.posi)
    dat<-afterposiExp
  }
  if(batch==TRUE){
            afterbatchExp<-ComBat(dat = dat, batch = batches, mod = mod, par.prior = par.prior, prior.plots = prior.plots, mean.only = mean.only.batch)
            dat<-afterbatchExp
  }
}







