plotDens1 <- function(mle, margin, int=NULL, col=1, lty=1, add=FALSE,
             xlim=NULL, ylim=NULL, xlab="", ylab="", main="", sub=""){

    p <- mle$p
    rects <- mle$rects

    # Check input
    if (!is.vector(p)) {
        stop("invalid argument 'mle$p'")
    }
    if (sum(p) > 1+1e-5 | sum(p) < 1-1e-5){
        warning("sum(mle$p) is not equal to 1")
    }
    if (!is.matrix(rects) | ncol(rects) != 4) {
        stop("invalid argument 'mle$rects': it must be a matrix with 4 columns")
    }
    if (sum(rects[, 1] > rects[, 2]) + sum(rects[, 3] > rects[, 4]) >= 1) {
        stop("invalid argument 'mle$rects': each row x1,x2,y1,y2 must satisfy x1<=x2 and y1<=y2")
    }
    if (length(p) != nrow(rects)) {
        stop("length(mle$p) must equal nrow(mle$rects)")
    }
    if (margin!=1 & margin!=2){
        stop("invalid argument 'margin'")
    }
    if (!is.null(int)){
        if (length(int)!=2 | (int[1]>=int[2]) ){
            stop("invalid argument 'int'")
        }
    }
    if (!is.logical(add)){
        stop("invalid argument 'add'")
    }

    m <- length(p) 
 
    # x-margin  
    if (margin==1){
       if (is.null(int)){
           intvl <- rects[,1:2]
       } else {
           ind.in <- c(1:m)[rects[,3]>=int[1] & rects[,4]<=int[2]]
           ind.hi <- c(1:m)[rects[,3]<int[1] & rects[,4]>int[1]] 
           ind.lo <- c(1:m)[rects[,4]>int[2] & rects[,3]<int[2]]
           intvl <- rects[c(ind.in, ind.hi, ind.lo),1:2]
           p <- c(p[ind.in], (int[2]-rects[,3])/(rects[,4]-rects[,3])*p[ind.hi], 
                  (rects[,4]-int[1])/(rects[,4]-rects[,3])*p[ind.lo])
       }
  
       if (sum(intvl[,1]==intvl[,2])>=1){
           stop("Cannot create density plot since some 
maximal intersections have equal x-coordinates") 
       }

       dens <- p / (intvl[,2]-intvl[,1])
       m <- length(p)       
       ord <- order(x <- c(intvl[,1],intvl[,2]), 
                    left <- c(rep(1,m),rep(0,m)), 
                    index <- c(c(1:m),c(1:m)))
       endpoints <- cbind(x[ord],left[ord],index[ord])

       x <- c(min(0,1000*min(mle$rects[,1])), endpoints[,1], 1000*max(mle$rects[,2]))       
       res <- c(rep(0,2*m+2))
       for (i in 1:(2*m)){
           r <- endpoints[i,3]
           if (endpoints[i,2]){
              res[i+1] <- res[i] + dens[r]
           } else { 
              res[i+1] <- res[i] - dens[r]
           }
       }
  
       if (is.null(xlim)){
           xlim <- range(rects[,1:2])
       }
       if (is.null(ylim)){
           ylim <- range(res)
       }
       if (!add){
           plot(x, res, type="s", col=col, lty=lty, xlim=xlim, ylim=ylim, 
                xlab=xlab, ylab=ylab, main=main, sub=sub)
       } else {
           lines(x, res, type="s", col=col, lty=lty)
       }
    }

    # y-margin  
    if (margin==2){
       if (is.null(int)){
           intvl <- rects[,3:4]
       } else {
           ind.in <- c(1:m)[rects[,1]>=int[1] & rects[,2]<=int[2]]
           ind.hi <- c(1:m)[rects[,1]<int[1] & rects[,2]>int[1]] 
           ind.lo <- c(1:m)[rects[,2]>int[2] & rects[,1]<int[2]]
           intvl <- rects[c(ind.in, ind.hi, ind.lo),3:4]
           p <- c(p[ind.in], (int[2]-rects[,1])/(rects[,2]-rects[,1])*p[ind.hi], 
                  (rects[,2]-int[1])/(rects[,2]-rects[,1])*p[ind.lo])
       }

       if (sum(intvl[,1]==intvl[,2])>=1){
           stop("Cannot create density plot since some 
maximal intersections have equal y-coordinates") 
       }

       dens <- p / (intvl[,2]-intvl[,1])
       m <- length(p)       
       ord <- order(x <- c(intvl[,1],intvl[,2]), 
                    left <- c(rep(1,m),rep(0,m)), 
                    index <- c(c(1:m),c(1:m)))
       endpoints <- cbind(x[ord],left[ord],index[ord])

       x <- c(min(0,1000*min(mle$rects[,1])), endpoints[,1], 1000*max(mle$rects[,2]))       
       res <- c(rep(0,2*m+2))
       for (i in 1:(2*m)){
           r <- endpoints[i,3]
           if (endpoints[i,2]){
              res[i+1] <- res[i] + dens[r]
           } else { 
              res[i+1] <- res[i] - dens[r]
           }
       }
  
       if (is.null(xlim)){
           xlim <- range(rects[,3:4])
       }
       if (is.null(ylim)){
           ylim <- range(res)
       }
       if (!add){
           plot(x, res, type="s", col=col, lty=lty, xlim=xlim, ylim=ylim, 
                xlab=xlab, ylab=ylab, main=main, sub=sub)
       } else {
           lines(x, res, type="s", col=col, lty=lty)
       }
   }
}  

