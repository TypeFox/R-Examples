#' Calculating effect size (Cohen's d) of independent two-sample t test
#'
#' @param m1          mean of sample 1
#' @param m2          mean of sample 2
#' @param sd1         standard deviation of sample 1
#' @param sd2         standard deviation of sample 2
#' @param n1          number of observations in sample 1
#' @param n2          number of observations in sample 2
#' @param t           t statistic
#' @param se1         standard error of sample 1
#' @param se2         standard error of sample 2
#' @param df          degree of freedom
#' @param alternative The test is two sided or one sided
#' @seealso     \code{\link{ES.t.one}}
#' @seealso     \code{\link{ES.t.paired}}
#' @export
#' @examples
#' ## mean, sd, n -> d
#' ES.t.two(m1=13.5,m2=5.5,sd1=4.1833,sd2=3.02765,n1=14,n2=10)
#'
#' ## mean se, n -> d
#' ES.t.two(m1=13.5,m2=5.5,se1=1.118034,se2=0.9574271,n1=14,n2=10)
#'
#' ## t and n -> d
#' ES.t.two(n1=14,n2=10,t=5.4349)
#' 
#' ## t, df and n -> d
#' ES.t.two(t = 5.4349, df = 21.982,n1=14,n2=10)
#' 
#' ## t and df -> d (assume n1=n2)
#' ES.t.two(t = 5.4349, df = 21.982)
ES.t.two <- function(m1=NULL,m2=NULL,sd1=NULL,sd2=NULL,n1=NULL,n2=NULL,t=NULL,se1=NULL,se2=NULL,df=NULL,alternative = c("two.sided", "one.sided")){
  d <- NULL
  alternative <- match.arg(alternative)
    if (sum(sapply(list(m1, m2, sd1, sd2, n1, n2), is.null)) == 0){
      sp <- sqrt(((n1-1)*sd1^2 + (n2-1)*sd2^2)/(n1+n2-2))
      d <- (m1-m2)/sp
    }else if(sum(sapply(list(n1, n2, t), is.null)) == 0){
      d <- t*sqrt((n1+n2)/n1/n2)
    }else if(sum(sapply(list(m1, m2, se1, se2, n1, n2), is.null)) == 0){
      sd1 <- se1 * sqrt(n1)
      sd2 <- se2 * sqrt(n2)
      sp <- sqrt(((n1-1)*sd1^2 + (n2-1)*sd2^2)/(n1+n2))
      d <- (m1-m2)/sp
    }else if(sum(sapply(list(t,df,n1,n2), is.null)) == 0){
      d <- (n1+n2)*t / sqrt(df*n1*n2)
    }else if(sum(sapply(list(t,df), is.null)) == 0){
      d <- 2*t / sqrt(df)
    }
  NOTE0="The alternative hypothesis is m1 > m2"
  if(alternative =="two.sided"){
    d <-abs(d)
    NOTE0="The alternative hypothesis is m1 != m2"
  }
  NOTE1="small effect size:  d = 0.2"
  NOTE2="medium effect size: d = 0.5"
  NOTE3="large effect size:  d = 0.8"
  NOTE=paste(NOTE0,NOTE1,NOTE2,NOTE3,sep="\n")
METHOD="effect size (Cohen's d) of independent two-sample t test"
structure(list(d = d, alternative = alternative, note=NOTE,method = METHOD), class = "power.htest")
}
