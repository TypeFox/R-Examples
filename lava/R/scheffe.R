##' Function to compute the Scheffe corrected confidence
##' interval for the regression line
##'
##' @title Calculate simultaneous confidence limits by Scheffe's method
##' @param model Linear model
##' @param newdata new data frame
##' @param conf.level confidence level (0.95)
##' @export
##' @examples
##' x <- rnorm(100)
##' d <- data.frame(y=rnorm(length(x),x),x=x)
##' l <- lm(y~x,d)
##' plot(y~x,d)
##' abline(l)
##' d0 <- data.frame(x=seq(-5,5,length.out=100))
##' d1 <- cbind(d0,predict(l,newdata=d0,interval="confidence"))
##' d2 <- cbind(d0,scheffe(l,d0))
##' lines(lwr~x,d1,lty=2,col="red")
##' lines(upr~x,d1,lty=2,col="red")
##' lines(lwr~x,d2,lty=2,col="blue")
##' lines(upr~x,d2,lty=2,col="blue")
scheffe <- function(model,newdata=model.frame(model),conf.level=0.95) {
    df <- model$df.residual
    p <- model$rank
    alpha <- 1-conf.level
    ## Scheffe value uses 1-tailed F critical value
    scheffe.crit <- sqrt(p*qf(1-alpha,p,df))
    ci <- predict(model,newdata,interval="confidence",level=conf.level)
    delta <- scheffe.crit/qt(1-alpha/2,df)
    ci[,2] <- ci[,1] -(ci[,1]-ci[,2])*delta
    ci[,3] <- ci[,1] +(ci[,3]-ci[,1])*delta
    return(ci)
}
