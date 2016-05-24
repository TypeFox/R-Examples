##spatial variation in the risk surface between different types
##return pointwise tolerance limits (PTLs) for plotting
mcseg.test <- function(pts, marks, h, stpts=NULL, ntest=100, proc=TRUE)
{
    adapt <- chkernel()
    ynames <- unique(marks)
    m <- length(ynames)
    ynames0 <- 1:m - 1
    names(ynames0) <- ynames
    y1 <- ynames0[as.character(marks)]
    nh <- length(h)
    stat<-rep(0, ntest)
    n <- length(y1)
    if(!is.null(stpts)) {
        nstpts <- nrow(stpts)
        tct <- matrix(NA, nrow=nstpts*m, ncol=ntest) ## m types????
    }
    c <- matrix(1, nrow=n, ncol=nh)
    alpha <- table(y1)/n
    for(i in 1:ntest){
        if(proc) cat("\rProcessing No.", i, "out of", ntest) 
        if(i==1) y2 <- y1 else y2 <- y1[sample(1:n)]
        if(nh > 1) {
            lcc <- .C("lcn", as.double(pts), as.integer(y2), as.integer(n),
                  as.double(h), as.integer(nh), as.integer(adapt$kernel),
                  as.double(c), lc=double(nh),
                  PACKAGE="spatialkernel")$lc
            ##ophndx <- which(lcc==max(lcc, na.rm=TRUE))
            ##oph <- h[ophndx[1]]
            oph <- h[which.max(lcc)]
        } else oph <- h
        p <- .C("hatpn", as.double(pts), as.integer(n), as.double(pts),
                as.integer(y2), as.integer(n), as.double(oph),
                as.integer(adapt$kernel),
                as.double(c), as.integer(m), p=double(n*m),
                PACKAGE="spatialkernel")$p
        p <- matrix(p, ncol=m)
        stat[i] <- sum(apply(p, 1, function(x) sum((x-alpha)^2)))
        if(!is.null(stpts)) {
            p <- .C("hatpn", as.double(stpts), as.integer(nstpts), as.double(pts),
                    as.integer(y2), as.integer(n), as.double(oph),
                    as.integer(adapt$kernel),
                    as.double(c), as.integer(m), p=double(nstpts*m),
                    PACKAGE="spatialkernel")$p
            tct[,i] <- (p-rep(alpha, each=nstpts))
        }
    }
    pv <- (ntest+1-rank(stat)[1])/ntest
    if(!is.null(stpts)) {
        pvtc <- apply(tct, 1, function(x) (ntest+1-rank(x)[1])/ntest)
        pvtc <- matrix(pvtc, ncol=m, dimnames = list(NULL, ynames))
        invisible(list(pvalue=pv, stpvalue=pvtc, pts=pts, marks=marks,
                       h=h, stpts=stpts, ntest=ntest))
    } else {
        invisible(list(pvalue=pv, pts=pts, marks=marks, h=h,
                       ntest=ntest))
    }
}
