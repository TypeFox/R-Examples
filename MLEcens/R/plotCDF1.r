plotCDF1 <- function(mle, margin, bound="b", int=NULL, surv=FALSE,
                     add=FALSE, col=1, lty=1, xlim=NULL, ylim=NULL, xlab="", 
                     ylab="", main="", sub=""){

    p <- mle$p
    rects <- mle$rects

    # Check input
    if (!is.vector(p)) 
        stop("invalid argument 'mle$p'")
    if (sum(p) > 1+1e-5 | sum(p) < 1-1e-5)
        warning("sum(mle$p) is not equal to 1")
    if (!is.matrix(rects) | ncol(rects) != 4) 
        stop("invalid argument 'mle$rects': it must be a 
matrix with 4 columns")
    if (sum(rects[, 1] > rects[, 2]) + sum(rects[, 3] > rects[, 4]) >= 1) {
        stop("invalid argument 'mle$rects': each row x1,x2,y1,y2 
must satisfy x1<=x2 and y1<=y2")
    }
    if (length(p) != nrow(rects)) 
        stop("length(mle$p) must equal nrow(mle$rects)")
    if (!is.logical(surv))
        stop("invalid argument 'surv'")
    if (margin!=1 & margin!=2)
        stop("invalid argument 'margin'")
    if (bound!="u" & bound!="l" & bound!="b")
        stop("invalid argument 'bound'")
    if (!is.logical(add))
        stop("invalid argument 'add'")
    if (!is.null(int)){
        if (length(int)!=2 | (int[1]>=int[2]) )
            stop("invalid argument 'int'")
    }

    m <- length(p)

    ## x-margin ##
    if (margin==1){
        if (is.null(xlim))
            xlim <- range(rects[,1:2])
        if ( (bound=="u") | (bound=="b") ){
            if (is.null(int)){
                ind <- c(1:m)
            } else {
                ind <- c(1:m)[rects[,3] > int[1] & rects[,3] <= int[2]]
            }
            x <- rects[ind,1]
            p <- mle$p[ind]
            ord <- order(x)
            x <- c(min(0, 1000*min(rects[,1])), x[ord],1000*max(rects[,2]))
            p <- c(0,p[ord],0)
            cp <- cumsum(p)
 
            if (add){
                if (surv){
                   lines(x, 1-cp, type="s", col=col, lty=lty)
                }else{ 
                   lines(x, cp, type="s", col=col, lty=lty)
                }
            } else {
                if (surv){
                   plot(x, 1-cp, type="s", col=col, lty=lty, xlim=xlim, 
                        ylim=ylim, xlab=xlab, ylab=ylab, main=main, sub=sub)
                }else{
                   plot(x, cp, type="s", col=col, lty=lty, xlim=xlim, 
                        ylim=ylim, xlab=xlab, ylab=ylab, main=main, sub=sub)
                }
            }
        }
        if ( (bound=="l") | (bound=="b") ){
            if (is.null(int)){
                ind <- c(1:m)
            } else {
                ind <- c(1:m)[rects[,4] > int[1] & rects[,4] <= int[2]]
            }
            x <- rects[ind,2]
            p <- mle$p[ind]
            ord <- order(x)
            x <- c(min(0,1000*min(rects[,1])), x[ord], 1000*max(rects[,2]))
            p <- c(0, p[ord],0)
            cp <- cumsum(p)
 
            if (add | (bound=="b") ){
                if (surv){
                   lines(x, 1-cp, type="s", col=col, lty=lty)
                }else{ 
                   lines(x, cp, type="s", col=col, lty=lty)
                }
            } else {
                if (surv){
                   plot(x, 1-cp, type="s", col=col, lty=lty, xlim=xlim, 
                        ylim=ylim, xlab=xlab, ylab=ylab, main=main, sub=sub)
                }else{
                   plot(x, cp, type="s", col=col, lty=lty, xlim=xlim, 
                        ylim=ylim, xlab=xlab, ylab=ylab, main=main, sub=sub)
                }
            }
        }
    }

    ## y-margin ##
    if (margin==2){
        if (is.null(xlim))
            xlim <- range(rects[,3:4])
 
        # upper bound:
        if ( (bound=="u") | (bound=="b") ){
            if (is.null(int)){
                ind <- c(1:m)
            } else {
                ind <- c(1:m)[rects[,1] > int[1] & rects[,1] <= int[2]]
            }
            y <- rects[ind,3]
            p <- mle$p[ind]
            ord <- order(y)
            y <- c(min(0, 1000*min(rects[,3])), y[ord], 1000*max(rects[,4]))
            p <- c(0,p[ord],0)
            cp <- cumsum(p)
 
            if (add){
                if (surv){
                   lines(y, 1-cp, type="s", col=col, lty=lty)
                }else{ 
                   lines(y, cp, type="s", col=col, lty=lty)
                }
            } else {
                if (surv){
                   plot(y, 1-cp, type="s", col=col, lty=lty, xlim=xlim, 
                        ylim=ylim, xlab=xlab, ylab=ylab, main=main, sub=sub)
                }else{
                   plot(y, cp, type="s", col=col, lty=lty, xlim=xlim, 
                        ylim=ylim, xlab=xlab, ylab=ylab, main=main, sub=sub)
                }
            }
        }
        if ( (bound=="l") | (bound=="b") ){
            if (is.null(int)){
                ind <- c(1:m)
            } else {
                ind <- c(1:m)[rects[,2] > int[1] & rects[,2] <= int[2]]
            }
            y <- rects[ind,4]
            p <- mle$p[ind]
            ord <- order(y)
            y <- c(min(0,1000*min(rects[,3])), y[ord], 1000*max(rects[,4]))
            p <- c(0, p[ord],0)
            cp <- cumsum(p)
 
            if (add | (bound=="b") ){
                if (surv){
                   lines(y, 1-cp, type="s", col=col, lty=lty)
                }else{ 
                   lines(y, cp, type="s", col=col, lty=lty)
                }
            } else {
                if (surv){
                   plot(y, 1-cp, type="s", col=col, lty=lty, xlim=xlim, 
                        ylim=ylim, xlab=xlab, ylab=ylab, main=main, sub=sub)
                }else{
                   plot(y, cp, type="s", col=col, lty=lty, xlim=xlim, 
                        ylim=ylim, xlab=xlab, ylab=ylab, main=main, sub=sub)
                }
            }
        }
    }
}
