#' Calculating effect size (Cohen's d) of paired two-sample t test
#'
#' @param md          mean difference (e.g., mean(x-y))
#' @param sd          standard deviation of mean differences (e.g., sd(x-y))
#' @param n           number of paires
#' @param t           t statistic
#' @param se          standard error of mean differences
#' @param df          degree of freedom
#' @param alternative The test is two sided or one sided
#' @seealso     \code{\link{ES.t.one}}
#' @seealso     \code{\link{ES.t.two}}
#' @export
#' @examples
#' ## md, sd -> d
#' ES.t.paired(md=-0.08062384,sd=1.401886)
#' 
#' ## md,se -> d
#' ES.t.paired(md=-0.08062384,se=0.1982566,n=50)
#' 
#' ## t, df -> d
#' ES.t.paired(t=-0.4067,df=49)
#' 
#' ## t, n -> d
#' ES.t.paired(t=-0.4067,n=50)
ES.t.paired <- function(md=NULL,sd=NULL,n=NULL,t=NULL,se=NULL,df=NULL,alternative = c("two.sided", "one.sided")){
  alternative <- match.arg(alternative)
  d <- NULL
    if(sum(sapply(list(md,sd), is.null)) == 0){
      d <- md/sd
    }else if(sum(sapply(list(md,se,n), is.null)) == 0){
      sd <- se * sqrt(n)
      d <- md/sd
    }else if(sum(sapply(list(t,df), is.null)) == 0){
      d <- t / sqrt(df)
    }else if(sum(sapply(list(t,n), is.null)) == 0){
      d <- t / sqrt(n-1)
    }
  NOTE0="The alternative hypothesis is md > 0"
  if(alternative =="two.sided"){
    d <-abs(d)
    NOTE0="The alternative hypothesis is md != 0"
  }
  NOTE1="small effect size:  d = 0.2"
  NOTE2="medium effect size: d = 0.5"
  NOTE3="large effect size:  d = 0.8"
  NOTE=paste(NOTE0,NOTE1,NOTE2,NOTE3,sep="\n")
  METHOD="effect size (Cohen's d) of paired two-sample t test"
  structure(list(d = d, alternative = alternative, note=NOTE,method = METHOD), class = "power.htest")
}
