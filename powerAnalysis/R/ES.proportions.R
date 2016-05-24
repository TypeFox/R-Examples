#' Compute effect size for a difference in proportions
#'
#' @param p1          Proportion of sample one
#' @param p2          Proportion of sample two or a constant proportion
#' @param alternative The test is two sided or one sided
#' @export
#' @examples
#' ES.proportions(0.65,0.45)
#' 
#' ES.proportions(0.25,0.05)
ES.proportions <- function(p1=NULL, p2=NULL,alternative = c("two.sided", "one.sided")){
  alternative <- match.arg(alternative)
  if(is.null(p1) || is.null(p2)){
    stop("'p1' and 'p2' can not be NULL")
  }
  if(p1 < 0 || p2 < 0 || p1>1 || p2>1){
    stop("'p0' and 'p1' must be in [0,1]")
  }
  h=2*(asin(sqrt(p1))-asin(sqrt(p2)))
  
  NOTE0="The alternative hypothesis is p1 > p2"
  if(alternative =="two.sided"){
    h <-abs(h)
    NOTE0="The alternative hypothesis is p1 != p2"
  }
  NOTE1="small effect size:  h = 0.2"
  NOTE2="medium effect size: h = 0.5"
  NOTE3="large effect size:  h = 0.8"
  NOTE=paste(NOTE0,NOTE1,NOTE2,NOTE3,sep="\n")
  METHOD="effect size of the difference in proportions"
  structure(list(h = h, p1 = p1, p2 = p2,note=NOTE,method=METHOD), class = "power.htest")                                     
}

