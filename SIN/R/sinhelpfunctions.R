#### HELP FUNCTIONS ##############
fisherz <- function(corrs){
  if(any(is.na(corrs))){
    return(NA)
  }
  if( (max(corrs)>1) || (min(corrs)< -1)){
    return("Argument is not a vector of correlations!")
  }
  else{
    return(0.5*log((1+corrs)/(1-corrs)))
  }
}
sdcor2cov <- function(stddev, corr){
    p <- (d <- dim(corr))[1]
    if (!is.numeric(corr) || length(d) != 2 || p != d[2]) 
        stop("`corr' is not a square numeric matrix")
    
    if (!is.numeric(stddev) || length(stddev) != d[2]) 
        stop("`stddev' and `corr' are not compatible")
    if (any(!is.finite(stddev))) 
        warning("stddev has non-finite entries")
    r <- corr
    r[] <- stddev * corr * rep(stddev, each = p)
    r
}
simpvalueVec <- function(corrs,n,p){
    temp <- sapply(corrs,fisherz)
    temp <- sapply(temp, abs)
    temp <- temp*sqrt(n)
    temp <- sapply(temp, pnorm)
    temp <- 1- ( 2*temp -1 )^(p*(p-1)/2)
    return(temp)
}
simpvalueMx <- function(corr,n,p){
  if(is.matrix(corr)){
    pp <- dim(corr)[1]
  }
  else{
    pp <- 1
  }
  if(pp==1){
    return(matrix(NA, 1,1))
  }
  else{
    temp <- simpvalueVec(c(corr),n,p)
    temp <- matrix(temp, pp,pp)
    diag(temp) <- NA
    return(temp)
  }
}
is.blocks <- function(blocks, p){
  if(!is.list(blocks)){
    return(FALSE)
  }
  if(!(all.equal(sort(unlist(blocks)), 1:p)==TRUE)){
    return(FALSE)
  }
  return(TRUE)
}
getgraph <- function(pvals, alpha, type="UG", blocks=NULL){
  getUG <- function(pvals, alpha){
    UG <- pvals
    diag(UG) <- 1
    UG <- matrix(as.numeric(UG <= alpha), ncol=dim(pvals)[1])
    dimnames(UG) <- dimnames(pvals)
    return(UG)
  }
  getDAG <- function(pvals, alpha){
    DAG <- getUG(pvals, alpha)
    DAG[lower.tri(DAG)] <- 0
    return(DAG)
  }
  getCG <- function(blocks, pvals, alpha){
    CG <- UG <- getUG(pvals,alpha)
    CG[lower.tri(CG)] <- 0
    for(i in 1:length(blocks)){
      CG[blocks[[i]],blocks[[i]]] <-
        UG[blocks[[i]],blocks[[i]]]
    }
    return(CG)
  }
  getBG <- function(pvals,alpha){
    return(2*getUG(pvals,alpha))
  }
  if( (!is.numeric(pvals)) || (!is.matrix(pvals)) ){
    return("pvals is not a matrix of p-values!")
  }
  if( (max(pvals[!is.na(pvals)])>1) ||
      (min(pvals[!is.na(pvals)])<0) ){
    return("pvals is not a matrix of p-values!")
  }
  if(!is.numeric(alpha)){
    return("alpha is not a significance level!")
  }
  if( (alpha<0) || (alpha>1) ){
    return("alpha is not a significance level!")
  }
  if(type=="CG"){
    if(!is.blocks(blocks, dim(pvals)[1])){
      return("blocks is not a valid block structure over the variables!")
    }
    return(getCG(blocks, pvals, alpha))
  }
  else{
    functioncall <- call(paste("get",type, sep=""), pvals, alpha)
    return(eval(functioncall))
  }
}
#### END FUNCTIONS ###########
