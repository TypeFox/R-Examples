#' Statistical Test with Wasserstein Metric
#' @param cases name of case group data (matrix sample * feature)
#' 
#' @param control names of control group data (matrix sample * feature)
#' 
#' @param test.stat test statistic
#' 
#' @param paranum the number of quatile discretization + 1. Default is discretized by 1 \%.
#' 
#' @param bsn the number of resampling. Default is bsn = 5000.
#' 
#' @param q power of Wasserstein metric. Default is q = 2.
#'
#' @param seed seed for random generator.
#' 
#' @examples
#' 
#' nrep <- 12
#' cases <- Map(rbeta,rep(30,nrep),rep(1,nrep),rep(5,nrep))
#' cases <- do.call("rbind",cases)
#' control <- Map(rbeta,rep(30,nrep),rep(1,nrep),rep(5,nrep))
#' control <- do.call("rbind",control)
#' d <- wasserMetric(cases,control)
#' testRes <- wasser.test(cases = cases,control = control,test.stat = d)

#' @author Yusuke Matsui & Teppei Shimamura
#' @references Yusuke Matsui, Masahiro Mizuta, Satoru Miyano and Teppei Shimamura.(2015) D3M:Detection of differential distributions of methylation patterns (submitted). BIORXIV/2015/023879.
#' @references Antonio Irpino and Rossanna Verde.(2015) Basic Statistics for distributional symbolic variables: a new metric-based approach. Adv.Data.Anal.Classif(9) 143--175
#' @return  list of p-value and test statistics.
#' 
#' @export

wasser.test <- function(cases, control,test.stat, paranum = 101, bsn = 5000, q = 2,seed = 100){
  #library(D3M)
  
  #library(Rcpp)
  
  set.seed(seed)
  
  nsample <- ncol(cases) + ncol(control)
  
  shuffleID <- sapply(1:bsn,function(j)sample(nsample,nsample,replace=FALSE))
  
  res <- permCpp(casesMat = cases, controlMat = control, shuffleID = shuffleID, bsn = bsn, qn = paranum, d = test.stat, q = q)

  o <- list(pval = res[[1]], tets.stat = test.stat)
  
  return(o)
}
