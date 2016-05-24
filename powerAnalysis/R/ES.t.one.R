#' Calculating effect size (Cohen's d) of one-sample t test
#'
#' @param m           mean of sample
#' @param sd          standard deviation of sample
#' @param n           number of observations
#' @param t           t statistic
#' @param se          standard error of sample 1
#' @param df          degree of freedom
#' @param mu          population mean
#' @param alternative The test is two sided or one sided
#' @seealso     \code{\link{ES.t.two}}
#' @seealso     \code{\link{ES.t.paired}}
#' @export
#' @examples
#' ## mean, sd and mu -> d
#' ES.t.one(m=-0.0938268,sd=0.9836668,mu=0)
#' 
#' ## mean, se, n and mu -> d
#' ES.t.one(m=-0.0938268,se=0.1391115,n=50,mu=0)
#' 
#' ## t and df -> d (df=n-1)
#' ES.t.one(t = -0.6745,df = 49)
#' 
#' ## t and n -> d ((df=n-1))
#' ES.t.one(t = -0.6745,n = 50)
ES.t.one <- function(m=NULL,sd=NULL,n=NULL,t=NULL,se=NULL,df=NULL,mu=NULL,alternative = c("two.sided", "one.sided")){
  alternative <- match.arg(alternative)
  d <- NULL
    if(sum(sapply(list(m,sd,mu), is.null)) == 0){
      d <- (m-mu)/sd
    }else if(sum(sapply(list(m,se,mu,n), is.null)) == 0){
      sd=se * sqrt(n)
      d <- (m-mu)/sd
    }else if(sum(sapply(list(t,df), is.null)) == 0){
      d <- t / sqrt(df)
    }else if(sum(sapply(list(t,n), is.null)) == 0){
      df=n-1
      d <- t / sqrt(df)
    }
  NOTE0="The alternative hypothesis is m > mu"
  if(alternative =="two.sided"){
    d <-abs(d)
    NOTE0="The alternative hypothesis is m != mu"
  }
  NOTE1="small effect size:  d = 0.2"
  NOTE2="medium effect size: d = 0.5"
  NOTE3="large effect size:  d = 0.8"
  NOTE=paste(NOTE0,NOTE1,NOTE2,NOTE3,sep="\n")
  METHOD="effect size (Cohen's d) of one-sample t test"
  structure(list(d = d, alternative = alternative, note=NOTE,method = METHOD), class = "power.htest")
}
