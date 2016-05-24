qqfun <- function(x, distribution="norm", ylab=deparse(substitute(x)),
            xlab=paste(distribution, "quantiles"), main=NULL, las=par("las"),
            envelope=.95, labels=FALSE, col=palette()[4], lcol=palette()[2], 
            xlim=NULL, ylim=NULL, lwd=1, pch=1, bg=palette()[4], cex=.4,
            line=c("quartiles", "robust", "none"), ...)
{
    result <- NULL
    line <- match.arg(line)
    good <- !is.na(x)
    ord <- order(x[good])
    ord.x <- x[good][ord]
    qfun <- eval(parse(text=paste("q",distribution,sep="")))
    dfun <- eval(parse(text=paste("d",distribution,sep="")))
    n <- length(ord.x)
    P <- ppoints(n)
    z <- qfun(P, ...)
    plot(z, ord.x, xlab=xlab, ylab=ylab, main=main, las=las, col=col, pch=pch, cex=cex, bg=bg, xlim=xlim, ylim=ylim)
    if (line=="quartiles") {
        Qx <- quantile(ord.x, c(.25,.75))
        Qz <- qfun(c(.25,.75), ...)
        b <- (Qx[2]-Qx[1])/(Qz[2]-Qz[1])
        a <- Qx[1]-b*Qz[1]
        abline(a, b, col=lcol, lwd=lwd)
    }
    if (line=="robust") {
        for(p in c("MASS")) {
           if (length(grep(paste("^package:", p, "$", sep=""), search())) == 0) {
              if (!require(p, quietly = TRUE, character.only=TRUE))
              warning(paste("qqfun needs package `", p, "' to be fully functional; please install", sep=""))
           }
        }
        coef <- coefficients(MASS::rlm(ord.x~z))
        a <- coef[1]
        b <- coef[2]
        abline(a,b,col=palette()[2])
    }
    if (line != 'none' & envelope != FALSE) {
        zz <- qnorm(1-(1-envelope)/2)
        SE <- (b/dfun(z, ...))*sqrt(P*(1-P)/n)
        fit.value <- a+b*z
        upper <- fit.value+zz*SE
        lower <- fit.value-zz*SE
        lines(z, upper, lty=2, lwd=lwd/2, col=lcol)
        lines(z, lower, lty=2, lwd=lwd/2, col=lcol)
    }
    if (labels[1]==TRUE & length(labels)==1) labels<-seq(along=z)
    if (labels[1] != FALSE) {
        selected <- identify(z, ord.x, labels[good][ord])
        result <- seq(along=x)[good][ord][selected]
    }
    if (is.null(result)) invisible(result) else sort(result)
}
