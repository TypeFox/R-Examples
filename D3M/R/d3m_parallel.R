#' Parallel Version of Two Sample Test with Distribution-Valued Data
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
#' @param n.core the number of cores to be used.
#' 
#' @param seed seed for random number generator assigned for each cluster. Default is c(100,200,300,400) assuming n.core=4.
#' 
#' @param type choose the cluster type. "SOCK","MPI","PVM"or "NWS". Default is "SOCK". More detail information, see ?clusterMap in "parllel" and/or "snow" package.
#' 
#' @return result of d3m in the list format, including p-value, test statistics, data of cases, data of control. Each element of list correponds to result from each core.
#' 
#' @author Yusuke Matsui & Teppei Shimamura
#' @references Yusuke Matsui, Masahiro Mizuta, Satoru Miyano and Teppei Shimamura.(2015) D3M:Detection of differential distributions of methylation patterns (submitted). BIORXIV/2015/023879.
#' @references Antonio Irpino and Rossanna Verde.(2015) Basic Statistics for distributional symbolic variables: a new metric-based approach. Adv.Data.Anal.Classif(9) 143--175
#' @examples
#' library(D3M)
#' nrep <- 12
#' cases <- Map(rbeta,rep(30,nrep),rep(1,nrep),rep(5,nrep))
#' cases <- do.call("rbind",cases)
#' control <- Map(rbeta,rep(30,nrep),rep(1,nrep),rep(5,nrep))
#' control <- do.call("rbind",control)
#' ## do not run.
#' #d3m.parallel(cases,control,rm.mean = FALSE, rm.var = FALSE, bsn = 1000)
#' 
#' @export
#' 


d3m.parallel <- function(cases, control, rm.mean = F, rm.var = F, paranum = 101 ,q = 2, bsn = 5000, n.core = 4, type = "SOCK", seed = c(100,200,300,400)){
  
  #library(snow)
  #library(Rcpp)
  
  
  cl <- parallel::makeCluster(n.core,type = type)
  
  worker.init <- function(packages){
    for(p in packages){
      library(p, character.only=TRUE)
    }
    NULL
  }
  
  #snow::clusterCall(cl,worker.init,packages)
  parallel::clusterCall(cl,worker.init, 'D3M')
  
  if(rm.mean & rm.var){
    
    cases <- scale(x = cases,center = T,scale = T)
    
    control <- scale(x = control,center = T,scale = T)
    
  }else if(rm.mean & !rm.var){
    
    control <- scale(x = control,center = T,scale = F)
    
  }else if(!rm.mean & rm.var){
    
    cases <- scale(x = cases,center = F,scale = T)
  }
  
  #nrep <- 1000
  #cases <- Map(rbeta,rep(100,nrep),rep(1,nrep),rep(5,nrep))
  #cases <- do.call("rbind",cases)
  #control <- Map(rbeta,rep(100,nrep),rep(1,nrep),rep(5,nrep))
  #control <- do.call("rbind",control)
  
  ncases<- ncol(cases)
  
  ncontrol<- ncol(cases)
  
  colnames(cases) <- as.character(1:ncol(cases))
  
  colnames(control) <- as.character(1:ncol(control))
  
  index <- seq(1,nrow(cases),1)
  
  cl_index <- parallel::clusterSplit(cl,index)
  
  vec <- vector("list",n.core)
  
  
  
  for(i in seq_along(cl_index)){
    
    vec[[i]] <- rep(i,length(cl_index[[i]]))
    
  }
  
  if(is.null(seed)){
    seed <- as.integer(seq(100, 1000, length = n.core))
  }
  
  vec <- unlist(vec)
  
  cases <- cbind(cases, cind = vec)
  
  control <- cbind(control, cind = vec)
  
  subcases <- lapply(1:n.core, function(i)subset(cases, cases[,"cind"] == i)[, -(ncases+1)])
  
  subcontrol <- lapply(1:n.core, function(i)subset(control, control[,"cind"]==i)[, -(ncontrol+1)])
  
  res <- parallel::clusterMap(cl = cl, fun = d3m, subcases, subcontrol, paranum, q, bsn, seed)
  
  parallel::stopCluster(cl)
  
  return(res)
}
