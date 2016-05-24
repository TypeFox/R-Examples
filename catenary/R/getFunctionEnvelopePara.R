#' Function to return function envelope for parabola
#' 
#' Use boostrap to get bands of possible fits to data using quadratic
#'
#' @param data data frame with columns \code{x} and \code{y}
#' @return data frame with x, lwr and upr
#' @author Jono Tuke, Matthew Roughan
#' @export
#' @keywords internal
#' @note February 12 2013
#' @examples
#' x <- runif(100,0,10)
#' y <- x^2 + 2*x + 3 + rnorm(100,sd=20)
#' df <- data.frame(x=x,y=y)
#' plot(y~x,data=df,pch=16,cex=0.5)
#' bounds <- getFunctionEnvelopePara(data=df,x=seq(0,10,l=100))
#' lines(bounds$x,bounds$lwr)
#' lines(bounds$x,bounds$upr)
getFunctionEnvelopePara <- function(data,R=1000,x){
  # Create boot stat function
  fn <- function(data,index,x){
    data <- data[index,]
    tmp.lm <- lm(y~ x + I(x^2),data=data)
    return(predict(tmp.lm,newdata=data.frame(x=x)))
  }
  # Fit boot strap
  b <- boot::boot(data=data,statistic=fn,R=R,x=x)
  # b contains all the sample predicted values so get this
  tmp <- apply(b$t,2,quantile,probs=c(0.025,0.975))
  bounds <- data.frame(x=x,lwr=tmp[1,],upr=tmp[2,])
  # return the df frame with lower and upper
  return(bounds)
}

