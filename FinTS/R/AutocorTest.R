AutocorTest <- function(x, lag=ceiling(log(length(x))),
            type=c("Ljung-Box", "Box-Pierce", "rank"),
            df=lag ){
  type <- match.arg(type)
  if(type=="rank"){
    x <- rank(x) 
    type <- "Ljung-Box"
  }
#
  LjB <- Box.test(x, lag, type)
#
  if(df<=0) df <- 1
  if(df != lag){
    LjB$parameter <- c(df=df) 
#   Correct the p.value  
    LjB$p.value <- pchisq(LjB$statistic, df, lower.tail=FALSE) 
    LjB$method <- paste(LjB$method, " (lag = ", lag, ")", sep="")
  }
  LjB$data.name <- deparse(substitute(x))
  LjB$Total.observ <- length(x)
#  
  LjB  
}
