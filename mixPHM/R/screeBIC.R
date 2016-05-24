`screeBIC` <-
function(x, lty = 1, col = NA, pch = 19, type = "b", 
main = "BIC Screeplot", xlab = "Number of Components", ylab = "BIC", legpos = "topright",...)
{
# x object of class BICmat
# produces a screeplot of BICmat
# additional graphical parameters can be assigned (...)

if (class(x)!="BICmat") {stop("Object of class BICmat (from function msBIC) needed!")}

lK <- max(x$K)
lx <- dim(x$BICmat)[1]
if (any(is.na(col))) col <- 1:lx

matplot(x$K,t(x$BICmat),type=type,lty=lty,col=col,main=main,xlab=xlab,
        ylab=ylab,pch=19,xaxp = c(min(x$K),max(x$K),length(x$K)-1),...)
#legend(legpos,lty=lty,col=col,legend=x$method)
}

