#' Calculating effect size (Cohen's f) of one-way anova for means with equal observations in each group 
#'
#' @param data a matrix or data frame
#' @param sst  total sum of squares
#' @param ssb  sum of squares between groups
#' @export
#' @examples
#' set.seed(9);x=rnorm(50);y=rnorm(50)
#' z=rnorm(50);d=data.frame(x,y,z)
#' ES.anova.oneway(data=d)
#' 
#' ES.anova.oneway(sst=50,ssb=1)
ES.anova.oneway <- function(data=NULL,sst=NULL,ssb=NULL){
  f=NULL
  if(!is.null(data) && any(is.data.frame(data),is.matrix(data))){
    n=nrow(data)
    mc=colMeans(data)
    m=mean(mc)
    
    sst=sum((data-m)^2)
    ssb=sum(50*(mc-m)^2)
    f=sqrt(ssb/(sst-ssb))
  }else if(!is.null(sst) && !is.null(ssb) && all(sst>0,ssb>0)){
    f=sqrt(ssb/(sst-ssb))
  }else{
    stop("'data', or both 'sst' and 'ssb' should be provided\n  'sst' and 'ssb' should be nonnegative\n  'data' should be a data frame or a matrix")
  }
  NOTE1="small effect size:  f = 0.1"
  NOTE2="      medium effect size: f = 0.25"
  NOTE3="      large effect size:  f = 0.4"
  NOTE=paste(NOTE1,NOTE2,NOTE3,sep="\n")
  METHOD="effect size (Cohen's f) of one-way anova test for means"
  structure(list(f = f, sst=sst,ssb=ssb,note=NOTE, method=METHOD), class = "power.htest")
}

