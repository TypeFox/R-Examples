"simulateLinear" <-
function(sd = 2, npoints=5, nrep=4, nsets=200, type="xy", seed=21)
{
    if(!is.null(seed))set.seed(seed)
    nval <- npoints*nrep
    tmp <- data.frame(x = rep(1:npoints, rep(nrep, npoints)))
    p.aov <- array(0, nsets)
    p.slope <- array(0, nsets)
        for(i in 1:nsets) {
        tmp$y <- 100 + 0.8 * tmp$x + rnorm(nval, 0, sd)
        u <- lm(y ~ factor(x), data = tmp)
        z <- summary.aov(u)
        p.aov[i] <- z[[1]][1,"Pr(>F)"]
        u <- lm(y ~ x, data = tmp)
        z1 <- summary(u)
        p.slope[i] <- z1$coef[2, 4]
     }
     logit <- function(p)log(p/(1-p))
     x <- logit(p.aov)
     y <- logit(p.slope)
     xlim <- range(c(x,y), na.rm = TRUE)
     if(type=="xy"){
        oldpar <- par(mar = par()$mar - c(.5, 0, 2, 0), mgp = c(2.75, 0.5, 0))
        on.exit(par(oldpar))
        plot(x, y, xlim=xlim, ylim=xlim, xlab="", ylab="", cex=0.75, 
            axes=FALSE, main="")
        pval <- c(0.001, 0.01, 0.1, 0.5, 0.9)
        xpos <- logit(pval)
        axis(1, at=xpos, labels=paste(pval))
        axis(2, at=xpos, labels=paste(pval))
        box()
        mtext(side=1, line=2.5, "p-value: Qualitative aov comparison")
        mtext(side=2, line=2.5, "p-value: Test for linear trend")    
        abline(0, 1)
        } else 
        if(type=="density"){
        oldpar <- par(mfrow=c(1,2), mar = par()$mar - c(.5, 0, 2, 0), mgp = c(2.75, 0.5, 0))
        on.exit(par(oldpar))
        denx <- density(x)
        deny <- density(y)
        ylim <- c(0, max(c(denx$y, deny$y)))
        plot(denx, type="l", xlim=xlim, ylim = ylim, axes=FALSE, yaxs="i", main="",
            xlab="Density curves - 2 sets of p-values")
        topleft <- par()$usr[c(1,4)]
        legend(x=topleft[1], y=topleft[2], lty=c(1,2), legend=c("aov","linear"), bty="n")
        pval <- c(0.001, 0.01, 0.1, 0.5, 0.9)
        xpos <- logit(pval)
        axis(1, at=xpos, labels=paste(pval))
        lines(deny, lty=2)
        plot(density(x-y), main="", xlab="Difference in p-values, logit scale", 
            bty="n", yaxs="i")
        axis(1)
        }
     frac <- sum(p.slope<p.aov)/nsets
     cat("\nProportion of datasets where linear p-value < aov p-value =", frac, "\n")
    invisible()
}
