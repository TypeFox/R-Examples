qqnormp<-function(y, ylim, p, main, xlab, ylab, ...) UseMethod("qqnormp")

qqnormp.default<-function(y, ylim, p=2, main="Exponential Power Distribution Q-Q plot",
            xlab="Theoretical Quantiles", ylab="Sample Quantiles",...) 
{
y<-y[!is.na(y)]
if (0 == (n<-length(y))) stop("y is empty")
if (missing(ylim)) ylim<-range(y)
x<-qnormp(ppoints(n),p=p)[order(order(y))]
plot(x,y,main=main, xlab=xlab, ylab=ylab, ylim=ylim, ...)
mtext(paste("p=",p), 3, 0.25)
invisible(list(x=x,y=y))
}

qqlinep<-function(y,p=2,...)
{
    y <- quantile(y[!is.na(y)],c(0.25, 0.75))
    x <- qnormp(c(0.25, 0.75),p=p)
    slope <- diff(y)/diff(x)
    int <- y[1]-slope*x[1]
    abline(int, slope, ...)
}

