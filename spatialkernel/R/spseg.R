##use data within a polygon 
##select bandwidth, calculate phat, and do MC test 
##opt=1,2,3 to do task one, two or three
spseg <- function(pts, marks, h, opt=2, ntest=100, poly=NULL, 
  delta=min(apply(apply(pts, 2, range), 2, diff))/100, proc=TRUE)
{
    if(!(opt %in% 1:3)) stop("\nargument opt must be one of 1, 2, or 3.\n")
    ans <- list(pts=pts, marks=marks, h=h, opt=opt)
    ## opt==1
    if(length(h) > 1) {
      if(proc) cat("\nCalculating cross-validated likelihood function\n")
      cv <- cvlogl(pts, marks, h)$cv
      hopt <- h[which.max(cv)]
      ans <- c(list(cv=cv, hcv=hopt), ans)
    } else {
      hopt <- h ## for opt=2, phat
    }
    if(opt >= 2) { ##phat and mcseg.test
      if(!is.null(poly)) ans$poly <- poly
      ## opt==2
      if(proc) cat("\nCalculating type-specific probabilities\n")
      mtypes <- unique(marks)
      m <- length(mtypes)
      if(is.null(poly)) {
        xyrng <- apply(pts, 2, range)
      } else xyrng <- apply(poly, 2, range)
      gridxy <- list(gridx=seq(xyrng[1,1]+delta/2, xyrng[2,1]-delta/2, by=delta),
        gridy=seq(xyrng[1,2]+delta/2, xyrng[2,2]-delta/2, by=delta))
      gridpts <- as.matrix(expand.grid(gridxy))
      ngrid <- nrow(gridpts)
      if(is.null(poly)) gridndx <- rep(TRUE, 1:ngrid) else 
        gridndx <- which(pinpoly(poly, gridpts)>0)
      p <- matrix(NA, ncol=m, nrow=ngrid)
      tmp <- phat(gridpts[gridndx,], pts, marks, hopt)$p
      p[gridndx,] <- tmp; colnames(p) <- colnames(tmp)
      ans <- c(gridxy, list(p=p), ans)
      if(opt==3) {
        ## opt==3
        if(proc) cat("\nMonte Carlo testing\n")
        stp <- matrix(NA, ncol=m, nrow=ngrid)
        mc <- mcseg.test(pts, marks, h, stpts=gridpts[gridndx,], 
          ntest=ntest, proc=proc)[c("pvalue", "stpvalue")]
        stp[gridndx,] <- mc$stpvalue; colnames(stp) <- colnames(mc$stpvalue)
        ans <- c(list(pvalue=mc$pvalue, stpvalue=stp, ntest=ntest), ans)
      }
    }
    if(proc) cat("\n")
    ans
}

plotcv <- function(obj, ...) plot(obj$h, obj$cv, type="l", ...)

plotphat <- function(obj, types=unique(obj$marks), sup=TRUE, col=risk.colors(10), 
    breaks=seq(0,1,length=length(col)+1), ...) 
{
    if(obj$opt<=1) stop("\nRun phatmctest() with argument opt=2 or 3.\n")
    sapply(types, function(x) match.arg(x, unique(obj$marks)))
    m <- length(types)
    for(j in 1:m) {
        if(is.null(obj$poly)) {
            plot(obj$pts[1,], xlab="", ylab="", xlim=range(obj$pts[,1]), ylim=range(obj$pts[,2]),
                 asp=1, main=types[j], type="n")
        } else {
            plot(obj$pts[1,], xlab="", ylab="", xlim=range(obj$poly[,1]), 
                 ylim=range(obj$poly[,2]), asp=1, main=types[j], type="n")
        }
        image(obj$gridx, obj$gridy, matrix(obj$p[,types[j]], ncol=length(obj$gridy)), 
              zlim=0:1, add=TRUE, col=col, breaks=breaks) 
        if(!is.null(obj$poly)) {
            polygon(obj$poly)
            if(sup) {
                ndx <- which(obj$marks == types[j]) ##not gridpts, pts!!!
                points(obj$pts[ndx,], ...)
            }
        }
        if(m>1 && j<m && interactive()) readline("\nPress Enter to plot next one ...")
    }
}

plotmc <- function(obj, types=unique(obj$marks), quan=c(0.05, 0.95), sup=FALSE, 
    col=risk.colors(10), breaks=seq(0,1,length=length(col)+1), ...) 
{
    if(obj$opt<=2) stop("\nRun phatmctest() with argument opt=3 first.\n")
    sapply(types, function(x) match.arg(x, unique(obj$marks)))
    m <- length(types)
    for(j in 1:m) {
        if(is.null(obj$poly)) {
            plot(obj$pts[1,], xlab="", ylab="", xlim=range(obj$pts[,1]), ylim=range(obj$pts[,2]),
                 asp=1, main=types[j], type="n")
        } else {
            plot(obj$pts[1,], xlab="", ylab="", xlim=range(obj$poly[,1]), 
                 ylim=range(obj$poly[,2]), asp=1, main=types[j], type="n")
        }
        image(obj$gridx, obj$gridy, matrix(obj$p[,types[j]], ncol=length(obj$gridy)), 
              zlim=0:1, add=TRUE, col=col, breaks=breaks)
        if(!is.null(quan)) {
            contour(obj$gridx, obj$gridy,
                    matrix(obj$stpvalue[,types[j]], ncol=length(obj$gridy)),
                    levels=quan, add=TRUE, ...)
        } 
        if(!is.null(obj$poly)) {
            polygon(obj$poly)
            if(sup) {
                ndx <- which(obj$marks == types[j]) ##not gridpts, pts!!!
                points(obj$pts[ndx,], ...)
            }
        }
        if(m>1 && j<m && interactive()) readline("\nPress Enter to plot next one ...")
    }
}
