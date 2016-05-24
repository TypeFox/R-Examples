#' Two Sample Test with Distribution-Valued Data
#' @param cases name of case group data (matrix)
#' 
#' @param control names of control group data (matrix)
#' 
#' @param rm.mean standarize each rows of cases and control to mean=0.
#' 
#' @param rm.var standarize each rows of cases and control to var=1.
#' 
#' @param paranum the number of quatile discretization + 1. Default is discretized by 1 \%.
#' 
#' @param q power of Wasserstein metric. Default is q = 2.
#' 
#' @param bsn the number of resampling. Default is bsn = 5000.
#' 
#' @param seed seed for random number generator.
#' 
#' @details this function is designed for two sample test based on Wasserstein metric. The function computes the the p-values based Wasserstein metric and resampling method. If rm.mean=F and rm.var=F, then statistical test is performed only based on more than 3rd order moments. 
#' 
#' @examples 
#' nrep <-12
#' cases <- Map(rbeta,rep(30,nrep),rep(1,nrep),rep(5,nrep)); cases <- do.call("rbind",cases)
#' control <- Map(rbeta,rep(30,nrep),rep(1,nrep),rep(5,nrep)); control <- do.call("rbind",control)
#' d3m(cases,control,paranum = 101, q = 2, bsn = 1000,seed = 100)
#' 
#' @return  pval p-value.
#' @return test.stat test statistic.
#' @return cases case group data used in the statistical test.
#' @return control control group data used in the statistical test.
#' @author Yusuke Matsui & Teppei Shimamura
#' @references Yusuke Matsui, Masahiro Mizuta, Satoru Miyano and Teppei Shimamura.(2015) D3M:Detection of differential distributions of methylation patterns (submitted). BIORXIV/2015/023879.
#' @references Antonio Irpino and Rossanna Verde.(2015) Basic Statistics for distributional symbolic variables: a new metric-based approach. Adv.Data.Anal.Classif(9) 143--175
#' @export
#' 

d3m <- function(cases, control, rm.mean = FALSE, rm.var = FALSE, paranum = 101, q = 2, bsn = 5000, seed = 100){
   
  #sourceCpp("./src/funcs.cpp")
  #library(Rcpp)
  
  if(rm.mean & rm.var){
    cases <- scale(x = cases,center = T,scale = T)
    control <- scale(x = control,center = T,scale = T)
  }else if(rm.mean & !rm.var){
    control <- scale(x = control,center = T,scale = F)
  }else if(!rm.mean & rm.var){
    cases <- scale(x = cases,center = F,scale = T)
  }

  d <- wasserMetric(cases,control,paranum = 101,q = 2)
  
  set.seed(seed)
  
  res <- wasser.test(cases,control,d,bsn = bsn)
  
  return(list(pval = res[[1]], test.stat = res[[2]], cases = cases, control = control))
  
}

