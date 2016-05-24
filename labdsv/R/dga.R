dga <- function(z,x,y,step=25,pres="+",abs="-",labcex=1,
    xlab=deparse(substitute(x)),ylab=deparse(substitute(y)),pch=1,title="",...)
{
    xstep <- seq(min(x),max(x),(max(x)-min(x))/step)
    ystep <- seq(min(y),max(y),(max(y)-min(y))/step)
    grid<-expand.grid(x=xstep,y=ystep)
    if (any(is.na(x))) {
        cat("Omitting plots with missing values \n")
        y <- y[!is.na(x)]
        z <- z[!is.na(x)]
        x <- x[!is.na(x)]
    }
    if (any(is.na(y))) {
        cat("Omitting plots with missing values \n")
        x <- y[!is.na(y)]
        z <- z[!is.na(y)]
        y <- x[!is.na(y)]
    }
    if (any(is.na(z))) {
        cat("Omitting plots with missing values \n")
        x <- y[!is.na(z)]
        y <- z[!is.na(z)]
        x <- x[!is.na(z)]
    }
    if (is.logical(z)) {
        cat(paste(" \n z = ",deparse(substitute(z)), 
            " \n x = ",deparse(substitute(x)), 
            " \n y = ",deparse(substitute(y)),"\n"))
        tmp.gam <- gam(z ~ s(x) + s(y),family=binomial)
        gam.pred <- matrix(predict.gam(tmp.gam,grid,type="response"),nrow=step+1)
        contour(xstep,ystep,gam.pred,levels=seq(0.2,0.8,0.2),labcex=labcex,
        xlab=xlab,ylab=ylab,main=title)
        points(x[z],y[z],pch=pres)
        points(x[!z],y[!z],pch=abs)
        attr(tmp.gam,'call') <- match.call()
        invisible(tmp.gam)
    } else {
        tmp.gam <- gam(z ~ s(x) + s(y),family=poisson)
        gam.pred <- matrix(predict.gam(tmp.gam,grid,type="response"),nrow=step+1)
        contour(xstep,ystep,gam.pred,labcex=1,
        xlab=xlab,ylab=ylab,main=title)
        quant <- quantile(z)
        points(x[z<=quant[2]],y[z<=quant[2]],cex=0.5,pch=pch)
        points(x[z>quant[2]&z<=quant[4]],
               y[z>quant[2]&z<=quant[4]],pch=pch)
        points(x[z>quant[4]],y[z>quant[4]],cex=1.5,pch=pch)
        attr(tmp.gam,'call') <- match.call()
        invisible(tmp.gam)
    }
}
