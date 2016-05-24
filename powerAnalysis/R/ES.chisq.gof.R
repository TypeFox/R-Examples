#' Compute effect size of chi-squared test of goodness of fit
#'
#' @param p1   a vector of frequencies or probabilities (alternative hypothesis). Frequencies will be rescaled to probabilities automatically. An error is given if any entry of p1 is negative.
#' @param p0   a vector of frequencies or probabilities of the same length of p1 (null hypothesis). Frequencies will be rescaled to probabilities automatically. An error is given if any entry of p0 is negative. Default value of p0 is a vector of 1/n with length of n. n is the length of p1.
#' @seealso    \code{\link{ES.chisq.assoc}}
#' @export
#' @examples
#' ES.chisq.gof(p1=c(10,20,30,40))
#' ES.chisq.gof(p1=c(0.1,0.2,0.3,0.4))
#' 
#' ES.chisq.gof(p1=c(10,20,30,40),p0=c(0.2,0.3,0.3,0.2))
#' ES.chisq.gof(p1=c(10,20,30,40),p0=c(20,30,30,20))
#' ES.chisq.gof(p1=c(0.1,0.2,0.3,0.4),p0=c(0.2,0.3,0.3,0.2))
#' ES.chisq.gof(p1=c(0.1,0.2,0.3,0.4),p0=c(20,30,30,20))
ES.chisq.gof <- function(p1=NULL, p0=rep(1/length(p1),length(p1))){
  if(is.null(p1)){
    stop("'p1' can not be NULL")
  }
  if(!is.vector(p1) || ! is.vector(p0)){
    stop("'p0' and 'p1' should be vector")
  }
  if(any(p1 < 0) || any(p0 < 0) || any(is.na(p1)) || any(is.na(p0))){
    stop("all entries of 'p0' and 'p1' must be nonnegative and finite")
  }
  if(length(p0) != length(p1)){
    stop("'p0' and 'p1' must have the same length")
  }
                                       
  if(all(p1>1)){
    p1=p1/sum(p1)
  }
  if(any(p1>1) || sum(p1)!=1){
    stop("probabilities in 'p1' must sum to 1")
  }
  if(all(p0>1)){
    p0=p0/sum(p0)
  }
  if(any(p0>1) || sum(p0)!=1){
    stop("probabilities in 'p0' must sum to 1")
  }
  w=sqrt(sum((p1-p0)^2/p0))
  
  if(w>0.5){
    effect="large"
  }else if(w>0.3){
    effect="medium"
  }else{
    effect="small"
  }
  NOTE0="Probabilities were same as provided or were rescaled from provided frequencies"
  NOTE1="small effect size:  w = 0.1"
  NOTE2="medium effect size: w = 0.3"
  NOTE3="large effect size:  w = 0.5"
  NOTE=paste(NOTE0,NOTE1,NOTE2,NOTE3,sep="\n")
  METHOD="effect size of chi-squared test of goodness of fit"
  structure(list(w = w, p1 = p1, p0 = p0,note=NOTE,method=METHOD), class = "power.htest")                                     
}

