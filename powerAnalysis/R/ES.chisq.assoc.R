#' Compute effect size of chi-squared test of association
#'
#' @param ct              a m x n Contingency Table (matrix with m rows and n columes)
#' @param chisq           the value the chi-squared test statistic
#' @param p               p value for the chi-squared test
#' @param df              degree of freedom (e.g., df=(m-1)*(n-1))
#' @param n               total number of observations (sample size)
#' @param mindf           the degrees of freedom for the variable with the smaller number of levels, if m > n, mindf=n-1, otherwise, mindf=m-1
#' @seealso               \code{\link{ES.chisq.gof}}
#' @export
#' @examples
#' counts <- matrix(c(225,125,85,95),nrow=2,byrow=TRUE)
#' ES.chisq.assoc(ct=counts)
#'
#' case <- c(225,85,100)
#' control <- c(125,95,125)
#' counts <- cbind(case,control)
#' ES.chisq.assoc(ct=counts)
#'
#' p1 <- c(225,85,100)
#' p2 <- c(125,95,125)
#' p3 <- c(175,90,113)
#' counts <- cbind(p1,p2,p3)
#' ES.chisq.assoc(ct=counts)
#'
#' ES.chisq.assoc(chisq=13.561,n=530,df=1,mindf=1)
#'
#' ES.chisq.assoc(p=0.000231,n=530,df=1,mindf=1)
ES.chisq.assoc <- function(ct=NULL, chisq=NULL, p=NULL, n=NULL, df=NULL, mindf=NULL){
  phi=NULL
  if(!is.null(ct)){
    row = nrow(ct)
    col = ncol(ct)
    if(any(ct < 0) || is.na(ct)){
      stop("all entries of 'ct' must be nonnegative and finite")
    }
    if(row < 2 || col <2){
      stop("'ct' must have at least two rows and two columes")
    }
    df = (row -1)*(col -1)
    n = sum(ct)
    mindf <- ifelse(row>col, col-1, row-1)
    chisq=chisq.test(ct)$statistic[[1]]
    p=chisq.test(ct)$p.value
    phi=sqrt(chisq/(n*mindf))
  }else if(!is.null(chisq)){
    if(chisq < 0){
      stop("chisq must be at least 0\n")
    }
    if(!is.null(n) && !is.null(df) && !is.null(mindf)){
      if(n<1){
        stop("total number of observations must be at least 1\n")
      }
      if(df < 1){
        stop("df must be at least 1\n")
      }
      if(mindf < 1){
        stop("mindf must be at least 1")
      }
      p=pchisq(chisq,df,lower.tail=FALSE)
      phi=sqrt(chisq/(n*mindf))
    }else{
      stop("n, df, and mindf are needed to calculate es(effect size)\n")
    }
  }else if(!is.null(p)){
    if(p < 0 || p > 1){
      stop("p must be in [0,1]\n")
    }
    if(!is.null(n) & !is.null(df) & !is.null(mindf)){
      if(n<1){
        stop("total number of observations must be at least 1\n")
      }
      if(df < 1){
        stop("df must be at least 1\n")
      }
      if(mindf < 1){
        stop("mindf must be at least 1\n")
      }
      chisq=qchisq(p,df,lower.tail=FALSE)
      phi=sqrt(chisq/(n*mindf))
    }else{
      stop("n, df, and mindf are needed to calculate es(effect size)\n")
    }
  }else{
    stop("one of ct, chisq and p is needed to calculate es(effect size)\n")
  } 
  METHOD="effect size of chi-squared test of association"
  
  large=0.5 * sqrt(mindf)
  medium=0.3 * sqrt(mindf)
  small=0.1 * sqrt(mindf)
  NOTE1=paste("small effect size:  phi = ",small,sep="")
  NOTE2=paste("      medium effect size: phi = ",medium,sep="")
  NOTE3=paste("      large effect size:  phi = ",large,sep="")
  NOTE=paste(NOTE1,NOTE2,NOTE3,sep="\n")
  structure(list(phi=phi, chisq=chisq, p=p, n=n, df=df, mindf=mindf,method=METHOD,note=NOTE), class = "power.htest") 
}
