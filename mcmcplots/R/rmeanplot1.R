rmeanplot1 <- function(x, col=NULL, lty=1, style=c("gray", "plain"), ...){
    x <- convert.mcmc.list(lapply(x, function(mco) apply(mco, 2, function(y) cumsum(y)/seq_along(y))))
    traplot1(x=x, col=col, lty=lty, style=style, xlab="Iteration", ylab="Running mean", ...)
}
