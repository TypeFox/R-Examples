#' Function to return function envelope for catenary
#' 
#' Use boostrap to get bands of possible fits to data using catenary
#'
#' @param data data frame with columns \code{x} and \code{y}
#' @param initial vector of starting values (c1,c2,lambda)
#' @return data frame with x, lwr and upr
#' @author Jono Tuke, Matthew Roughan
#' @export
#' @keywords internal
#' @note February 12 2013
#' @examples
#' x <- runif(100,-2,2)
#' y <- f(x=x,c1=1,c2=2,lambda=3) + rnorm(100)
#' df <- data.frame(x=x,y=y)
#' plot(y~x,data=df,pch=16,cex=0.5)
#' bounds <- getFunctionEnvelopeCat(data=df,initial=c(1,2,3),x=seq(-2,2,l=100))
#' lines(bounds$x,bounds$lwr)
#' lines(bounds$x,bounds$upr)
getFunctionEnvelopeCat <- function(data,R=1000,initial,x){
  # Create boot stat function
  fn <- function(data,index,x,initial){
    data <- data[index,]
    c1 <- initial[1];c2 <- initial[2]; lambda <- initial[3]
    pred <- rep(NA,length(x))
    tryCatch({obs.cat <- nls(y ~ c1 * cosh( (x-c2)/c1) + lambda,
                             data=data,
                             start=list(c1=c1,c2=c2,lambda=lambda))
              pred <- predict(obs.cat,newdata=data.frame(x=x))
    }, error = function(e){
      cat("cannot fit catenary\n")
    })
    return(pred)
  }
  # Fit boot strap
  b <- boot::boot(data=data,statistic=fn,R=R,x=x,initial=initial)
  # b contains all the sample predicted values so get this
  tmp <- apply(b$t,2,quantile,probs=c(0.025,0.975),na.rm=TRUE)
  bounds <- data.frame(x=x,lwr=tmp[1,],upr=tmp[2,])
  # return the df frame with lower and upper
  return(bounds)
}
