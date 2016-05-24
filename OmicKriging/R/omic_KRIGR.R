#' Run omic kriging on a set of correlation matrices and a given phenotype.
#' 
#' Universal kriging formula:
#'   lambda' = ( c + X m )' iSig
#'   m' = ( x - X' iSig c )' ( X' iSig X )^-1
#'   m' = ( t(x) - c' iSig X ) ( X' iSig X )^-1
#'   lambda' = (c' + m' X) iSig
#'   x: #covariates x ntest
#'   X: ntrain x #cov
#'   c: ntrain x ntest
#' 
#' @param corlist A list of correlation matrices used in Kriging. rownames and colnames
#'   of cor should be IID list and include idtest and idtrain.
#' @param H2vec has weights for each RM relatednes matrix
#' @param idtest A vector of sample IDs which constitute the test set.
#' @param idtrain A vector of sample IDs which constitute the training set.
#' @param pheno A data frame with rownames set as sample IDs and a column containing phenotype data.
#' @param phenoname The name of the column in pheno which contains phenotype data to test.
#' @param Xcova Data frame of covariates with rownames() set to sample IDs. 
#' 
#' @return A dataframe with three columns: sample ID, observed phenotype Ytest, and predicted phenotype Ypred 
#' 
#' @keywords prediction
#' 
#' @references Cressie 1993 Statistics for Spatial Data p.154
#'
#' @export
okriging <- function(idtest,idtrain=NULL,corlist,H2vec,pheno,phenoname,Xcova=NULL){
  idtest <- as.character(idtest)
  idtrain <- as.character(idtrain)
  nt <- length(idtest)
  nT <- length(idtrain)
  indall <- c(idtrain,idtest)
  if(length(unique(idtest))!=nt) warning('repeated test ids')
  if(length(unique(idtrain))!=nT) warning('repeated train ids')
  if(length(intersect(idtest,idtrain)>0)) warning('test id in training set')
  if(sum(H2vec<0) | sum(H2vec)>1) stop(' sum of weights > 1 or negative weights ')
  
  ## compute correlation matrix
  if(length(corlist)!=length(H2vec)) stop('number of correlation components (length(H2vec)) != number of corlist ')
  id <- diag(rep(1,nt+nT)) ## identity matrix
  Sigmall <- id * (1 - sum(H2vec))
  for(cc in 1:length(corlist)) Sigmall = Sigmall + H2vec[cc] * corlist[[cc]][indall,indall]

  ## row and colnames of cor should be IID
  if(sum(c(idtest,idtrain) %in% rownames(Sigmall))<(nt+nT)) stop('some correlations are missing')
  
  ## if no covariates, use intercept
  Xtest <- matrix(1,1,nt)
  Xtrain <- matrix(1,nT,1)

  if(!is.null(Xcova)) 
  {
    Xtest <- rbind(Xtest,matrix(t(Xcova[idtest,]),ncol(Xcova),nt))
    Xtrain <- cbind(Xtrain,as.matrix(Xcova[idtrain,]))
  }
  
  Ytrain <- pheno[idtrain,phenoname]
  
  ## iSig
  iSig <- solve( Sigmall[idtrain,idtrain] ) 
  
  ctvec <- matrix(Sigmall[idtest,idtrain],nt,nT )  ## correlation between new id and old id (nT x nt)
  cvec <- t(ctvec)
  mtvec <- (t(Xtest) - ctvec %*% iSig %*% Xtrain) %*% solve( t(Xtrain) %*% iSig %*% Xtrain)

  lambt <- (ctvec + mtvec %*% t(Xtrain)) %*% iSig
  Ypred <- lambt %*% Ytrain
  Ytest <- pheno[idtest,phenoname]
  res <- data.frame(IID=idtest,Ypred,Ytest)
  rownames(res) <- idtest
  return(res)
}

#' Multithreaded cross validation routine for Omic Kriging.
#'
#' This is a flexible cross validation routine which wraps the Omic Kriging
#' calculation. The user can specify the size of the test set, all the way to
#' "Leave One Out" cross validation. Additionally, all relevant  parameters in the
#' \code{\link{okriging}} function are exposed. This function uses the doParallel
#' package to distribute computation over multiple cores. If the phenotype is 
#' case/control, a ROCR AUC and GLM analysis is run and the results printed to screen.
#'
#' @param cor.list A list of correlation matrices used in Kriging. rownames and colnames
#'   of cor should be IID list and include idtest and idtrain.
#' @param h2.vec has weights for each RM relatednes matrix
#' @param pheno.df A data frame with rownames set as sample IDs and a column containing phenotype data.
#' @param pheno.id The name of the column in pheno which contains phenotype data to test.
#' @param covar.mat Data frame of covariates with rownames() set to sample IDs. 
#' @param nfold Select the number of cross validation rounds to run. The value "LOOCV"
#'   will run one round of cross validation for each sample in your dataset.
#'   The value "ncore" will set the test set size such that a single round
#'   runs on each core specified in the ncore option. Any numeric value
#'   will be set to the test size. Default runs 10 rounds of cross validation.
#' @param ncore The number of cores available to distribute computaition across
#'    If a numeric value is supplied, that number of cores is registered. If the
#'    value "all" is supplied, all available cores are used. 
#' @param verbose Report rounds on cross validation on standard out. 
#' @param ... Optional and unnamed arguments.
#'
#' @return  A dataframe with three columns: sample ID, observed phenotype Ytest, and predicted phenotype Ypred
#'
#' @keywords prediction, cross validation
#'
#' @include omic_KRIGR.R
#'
#' @import doParallel
#' @import parallel
#' @import foreach
#' @import ROCR
#' @importFrom utils flush.console
#' @importFrom stats glm 
#' @importFrom stats binomial 
#' @importFrom stats lm 
#' @export
krigr_cross_validation <- function(cor.list, pheno.df, pheno.id = 1, h2.vec, covar.mat = NULL, nfold = 10, ncore = "all", verbose = FALSE, ...) {
  ## dependencies
  ## functions
  '%&%' <- function(a, b) paste(a, b, sep="")

  ## split into groups based on the number of cores available
  rownames(pheno.df) <- pheno.df$IID
  sample.ids <- pheno.df$IID
  n.samples <- length(sample.ids)
  cat("Detected", n.samples, "samples...", "\n")    
  flush.console()
  
  ## detect cores  
  if(ncore == "all") {
    ncore <- detectCores()
  } else if (!is.numeric(ncore)) {
    stop("ncore supplied must either be numeric or the string \"all\".")
  }
  
  ## register cluster
  clust <- makeCluster(ncore)
  registerDoParallel(clust)

  ## set n-fold  
  if(nfold == "LOOCV") {
    nfold <- n.samples
    } else if(is.numeric(nfold)) {
    nfold <- nfold
    } else if(nfold == "ncore") {
    nfold <- ncore
    } else {
    nfold <- 10
    }
   
  ## print core and fold numbers
  if(nfold == n.samples) {
    cat("Set leave-one-out cross-validation...", "\n")
    flush.console()
    } else {
    cat("Set", nfold %&% "x", "cross-validation...", "\n")
    flush.console()
    }
  
  cat("With", ncore, "logical core(s)...", "\n") 
  flush.console()

  ## create groups
  if(nfold == n.samples) {
    rand.groups <- 1:n.samples
    } else {
    rand.groups <- sample(1:nfold, n.samples, replace = T)
    }
  group.df <- data.frame(rand.groups, sample.ids)
  colnames(group.df) <- c("group.id", "sample.id") 
 
  cat("Running OmicKriging...", "\n")
  flush.console()
  
  ## running kriging routine on each core for each testing group
  n.par <- unique(rand.groups)
  i <- 0 ## added i to the functions namespace so that R CMD Check catches it
  time <- system.time(
    res <- foreach(i = 1:length(n.par), .combine = rbind, .export = ls(envir=globalenv())) %dopar% {
      if (verbose) cat(n.par[i], "\n")
      flush.console()
    
      ## separate test/train for round i of the cross validation
      if(length(n.par) == n.samples){
        idtest <- as.character(sample.ids[i])
        } else {
        idtest <- as.character(group.df$sample.id[group.df$group.id == n.par[i]])
        }
      idtrain <- as.factor(group.df$sample.id[!(group.df$sample.id %in% idtest)])
      
      ## run kriging
      okriging(idtest, idtrain, 
            corlist <- cor.list, 
            H2vec <- h2.vec,
            pheno <- pheno.df,
            phenoname <- colnames(pheno.df)[pheno.id + 2],
            Xcova <- covar.mat
            )
      }
    )

  gc()

  ## summary
  if(length(unique(res$Ytest)) == 2) {
    auc <- function(predtype, phenotype){
      pred <- prediction(predtype, phenotype)
      perf <- performance(pred, "auc")
      aucval <- perf@y.values
      return(aucval)
      }
    cat("Summary of binary phenotype...", "\n")
    cat("Area under the ROC curve:", auc(res$Ypred,res$Ytest) %&% "...", "\n")
    convertpheno <- res$Ytest
    uniconvertpheno <- unique(convertpheno)
    convertpheno[convertpheno == uniconvertpheno[2]] = 0
    convertpheno[convertpheno == uniconvertpheno[1]] = 1
    sum <- summary(glm(convertpheno ~ res$Ypred, family = binomial))
    print(sum)
    } else {
    sum <- summary(lm(Ytest ~ Ypred, data = res))
    print(sum)
    }

  cat("Finished OmicKriging in", time[3], "seconds", "\n")
  return(res)

}
