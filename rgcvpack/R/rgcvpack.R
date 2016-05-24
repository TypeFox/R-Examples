#rgcvpack main interfaces

fitTps <- function(x, y, m=2, knots=NULL, scale.type="range", method="v",
                   lambda=NULL, cost=1, nstep.cv=80, verbose=FALSE, tau=0)
{
    this.call <- match.call()
    
    if (is.data.frame(x))  x <- as.matrix(x)
    if (is.data.frame(y))  y <- as.matrix(y)
    if (!is.numeric(x) || !is.numeric(y))
        stop("the design points x and observation y should be numeric")
    
    n <- length(y)
    if (n <= 0)
        stop("the observation y is not in the right form")
    
    if (length(x)==n) {#one dimensional case
        d <- 1
    } else {
        if (length(x)<=0 || length(x)%%n!=0)
            stop("the dim of x is not compatible with y")
        else {
            xdim <- dim(x)
            if (is.null(xdim) || length(xdim)!=2)
                stop("the design points x must be a 2 dimensional array")
            
            if (xdim[1]!=n)
                stop("the number of rows in x must be compatible with y")
            else
                d <- xdim[2]
        }
    }
    
    if (!is.numeric(m) || (m-round(m))!=0 || 2*m-d<=0)
        stop("the order m of the spline is incorrect")
    
    if (is.null(knots)) {
        full.knots <- TRUE
        knots <- x
    } else {
        if (is.data.frame(knots))  knots <- as.matrix(knots)
        
        if (d==1) {
            if (!is.numeric(knots) || length(knots)>n)
                stop("the knots specification is incorrect")
            else
                full.knots <- FALSE
        } else {
            if (!is.numeric(knots) || dim(knots)[1]>n || dim(knots)[2]!=d)
                stop("the knots specification is incorrect")
            else
                full.knots <- FALSE
        }
    }
    
    if (scale.type=="range") {
        xlbs <- apply(x, 2, min)
        xubs <- apply(x, 2, max)
        
        xs   <- aperm(apply(x, 1, function(x, lbs=xlbs, ubs=xubs)
                        { (x - lbs)/(ubs - lbs) }), c(2,1))
        ks   <- aperm(apply(knots, 1, function(x, lbs=xlbs, ubs=xubs)
                        { (x - lbs)/(ubs - lbs) }), c(2,1))
    } else if (scale.type=="none") {
        xs <- x
        ks <- knots
    } else {
        stop("the scale.type can only take 'range' or 'none'")
    }
    
    job <- ifelse(full.knots && any(duplicated(xs)), 1, 0)
    
    if (!is.character(method) || nchar(method)!=1) {
        stop("the method argument should be a character of length 1")
    } else {
        method <- tolower(method)
        if (method=="d") {
            if (is.null(lambda) || !is.numeric(lambda) || lambda<0)
                stop("the lambda argument is incorrect")
            
            lamlim <- log10(lambda)*c(1,1)
            job <- job + ifelse(full.knots, 100, 10)
        } else {
            if (method=="v") {
                lamlim <- numeric(2)
            }
            else
                stop("the method argument can only take 2 values (d or v)")
        }
    }
    
    if (!is.numeric(cost) || cost<=0)
        stop("the cost argument should be a postive number")

    if (!is.numeric(nstep.cv) || nstep.cv<=0)
        stop("the nstep.cv argument should be a postive number")
    
    if (!is.logical(verbose))
        stop("the argument verbose can only be TRUE/FALSE")
    
    if (!is.numeric(tau) || tau<0)
        stop("the argument tau should be nonnegative")
    
    job <- job + ifelse(!full.knots && tau>0, 1000, 0)
    
    if (verbose) {
        cat("x =\n")
        print(x)
        cat("y =\n")
        print(y)
        cat("n =", n, "\n")
        cat("m =", m, "\n")
        cat("d =", d, "\n")
        cat("knots =\n")
        print(knots)
        cat("method =", method, "\n")
        cat("lambda =", lambda, "\n")
        cat("nstep.cv =", nstep.cv, "\n")
        cat("tau =", tau, "\n")
        cat("job =", job, "\n")
    }
    
    #call the dtpss or dsnsm driver for tps fitting
    nnull <- choose(m+d-1, d)
    if (full.knots) {

        des <- xs; ybak <- y
        
        if (verbose) cat("\nfitting thin plate spline...")
        tpss.lst <- .Fortran("dtpss", des=as.double(des), lddes=as.integer(n),
                             nobs=as.integer(n), dim=as.integer(d),
                             m=as.integer(m), s=double(n), lds=as.integer(n),
                             ncov=as.integer(0), y=as.double(ybak),
                             ntbl=as.integer(nstep.cv), adiag=double(n),
                             lamlim=as.double(lamlim), cost=as.double(cost),
                             dout=double(5), iout=integer(4),
                             coef=double(nnull+n), svals=double(n),
                             tbl=double(nstep.cv*3),ldtbl=as.integer(nstep.cv),
                             auxtbl=double(3*3), work=double(n*(3+nnull+n)),
                             lwa=as.integer(n*(3+nnull+n)), iwork=integer(3*n),
                             liwa=as.integer(3*n), job=as.integer(job),
                             info=integer(1), DUP=FALSE,
                             PACKAGE="rgcvpack")[-c(21,22)]
        if (verbose) cat("Done.\n\n")
        
    } else {

        b <- nrow(knots); ybak <- y
        
        if (verbose) cat("fitting the thin plate spline...")        
        tpss.lst <- .Fortran("dtpkm", des=as.double(xs), lddes=as.integer(n),
                             nobs=as.integer(n), dim=as.integer(d),
                             m=as.integer(m), s=double(n), lds=as.integer(n),
                             ncov=as.integer(0), knots=as.double(ks),
                             ldknt=as.integer(b), nknt=as.integer(b),
                             y=as.double(ybak), ntbl=as.integer(nstep.cv),
                             adiag=double(n), lamlim=as.double(lamlim),
                             cost=as.double(cost), dout=double(5),
                             iout=integer(3), coef=double(nnull+b),
                             svals=double(b), tbl=double(nstep.cv*3),
                             ldtbl=as.integer(nstep.cv), auxtbl=double(3*3),
                             work=double(2*b*(1+b+n)+n),
                             lwa=as.integer(2*b*(1+b+n)+n), iwork=integer(2*n),
                             liwa=as.integer(2*n), tau=as.double(tau),
                             job=as.integer(job),
                             info=integer(1), DUP=FALSE,
                             PACKAGE="rgcvpack")[-c(24,25)]
        
        if (verbose) cat("Done.\n\n")
        
    }
    
    if (tpss.lst$info > 0)
        stop(paste("tps fitting exit with an error code", tpss.lst$info))
    
    coefd <- tpss.lst$coef[1:nnull]
    coefc <- tpss.lst$coef[-(1:nnull)]
    
    lamhat <- n*tpss.lst$dout[1]
    df <- (n - tpss.lst$dout[4])/cost
    gcv <- tpss.lst$auxtbl[4]
    gcv.grid <- as.data.frame(matrix(tpss.lst$tbl, nstep.cv, 3))
    names(gcv.grid) <- c("loglam", "fGCV", "PMSE")
    gcv.grid$fGCV[gcv.grid$fGCV==1e+20] <- NA
    object <- list(x=x, y=y, m=m, knots=knots, scale.type=scale.type,
                   method=method, lambda=lamhat, cost=cost, nstep.cv=nstep.cv,
                   tau=tau, df=df, gcv=gcv, xs=xs, ks=ks, c=coefc, d=coefd,
                   yhat=tpss.lst$y, svals=tpss.lst$svals, gcv.grid=gcv.grid,
                   call=this.call)
    class(object) <- "Tps"
    object
    
}

print.Tps <- function(x, digits=4, ...) {
    
    if (!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl)
    }
    
    m <- x$m; d <- ncol(x$x); res.df <- nrow(x$x) - x$df
    cat("\n Number of Observations:     ", nrow(x$x))
    cat("\n Null space dimension:       ", choose(m+d-1,d))
    cat("\n Model degrees of freedom:   ", format(signif(x$df,digits)))
    cat("\n Residual degrees of freedom:", format(signif(res.df,digits)))
    cat("\n GCV estimate of lambda:     ", format(signif(x$lambda)))
    cat("\n GCV score of the model:     ", format(signif(x$gcv)))
    cat("\n")
    invisible(x)
}

predict.Tps <- function(object, newdata=NULL, ...) {
    
    if (class(object)!="Tps") {
        stop("the predict.Tps function requires a Tps object")
    }
    
    if (is.null(newdata)) {
        return (object$yhat)
    }
    
    if (is.data.frame(newdata)) {
        newdata <- as.matrix(newdata)
    }
    
    ks <- object$ks
    d  <- ncol(ks)
    
    if (ncol(newdata)!=d) {
        stop("the newdata should have the same # of columns as basis")
    }

    if (object$scale.type=="range") {
        xlbs <- apply(object$x, 2, min)
        xubs <- apply(object$x, 2, max)
        
        newdats <- aperm(apply(newdata, 1, function(x, lbs=xlbs, ubs=xubs)
                           { (x - lbs)/(ubs - lbs) }), c(2,1))
    } else { #object$scale.type=="none"
        newdats <- newdata
    }
    
    npred  <- nrow(newdats)
    nknots <- nrow(ks)
    m      <- object$m
    nnull  <- choose(m+d-1, d)
    npar   <- nknots + nnull
    coef   <- c(object$d, object$c)
    lwa    <- npred*(nnull + nknots)
    
    pred <- .Fortran("dpred", pdes=as.double(newdats), ldpdes=as.integer(npred),
                     npred=as.integer(npred), dim=as.integer(d),m=as.integer(m),
                     desb=as.double(ks), lddesb=as.integer(nknots),
                     ndesb=as.integer(nknots), ps=double(npred),
                     ldps=as.integer(npred), ncov1=as.integer(0),
                     ncov2=as.integer(0), coef=as.double(coef),
                     npar=as.integer(npar), pred=double(npred),
                     work=double(lwa), lwa=as.integer(lwa),
                     iwork=integer(d), info=integer(1), DUP=FALSE,
                     PACKAGE="rgcvpack")$pred
    
    return(pred)
    
}

