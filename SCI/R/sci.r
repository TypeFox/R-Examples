.mon.vect <- function(first.mon,n){
### produce replicated sequences from 1:12 with the condition
### that the first value has the value of 'first.mon'
###
### first.mon : value in 1:12
### n : length of the resulting series
    if(!(first.mon%in%(1:12)))
        stop("'first.mon' should be in 1:12")
    mm <- (1:12)+(first.mon-1)
    mm[mm>12] <- mm[mm>12] - 12
    rep(mm,length.out=n)
}

fitSCI <- function(x,...)
    UseMethod("fitSCI")

fitSCI.default <- function(x,
                           first.mon,
                           time.scale,
                           distr,
                           p0,
                           p0.center.mass=FALSE, ## use Weibull plotting position function for p0 estimation
                           scaling=c("no","max","sd"),
                           mledist.par=list(),
                           start.fun=dist.start,
                           start.fun.fix=FALSE,
                           warn=TRUE,...){
    if(time.scale < 1)
        stop("time.scale must be > 1")
    if(time.scale%%1!=0)
        stop("time.scale must be a integer number (1, 2, ...)")
    x <- as.numeric(x)   
    scaling <- match.arg(scaling)
    scale.val <- switch(scaling,
                        no=1,
                        max=max(x,na.rm=TRUE),
                        sd=sd(x,na.rm=TRUE))
    x <- x/scale.val   
    nn <- length(x)
    time.index <- 1:nn
    mmon <- .mon.vect(first.mon=first.mon,n=nn)
    ## 1. backward moving average
    if(time.scale>1){
        ffilt <- rep(1,time.scale)/time.scale
        x <- as.numeric(filter(x,filter=ffilt,method="convolution",sides=1))
    }
    ## 2. find distribution parameter for each month
    x.fit <- vector("list", 12)
    names(x.fit) <- paste("M",1:12,sep="")
    x.fit.monitor <- rep(0,12)
    names(x.fit.monitor) <- names(x.fit)
    empty.fit <- try(start.fun(NA,distr),silent=TRUE)
    if(class(empty.fit)=="try-error"){
        stop("distribution ", distr," not defined for start.fun=",deparse(substitute(start.fun)),
             "\n or startfun not implemented correctly")
    } else {
        empty.fit <- unlist(empty.fit)
    }
    if(p0){
        empty.fit <- c(empty.fit,P0=NA)
    }
    mledist.par$distr <- distr
    for(mm in 1:12){
        mledist.par$data <- x[mmon==mm] ## select month
        mledist.par$data <- mledist.par$data[is.finite(mledist.par$data)]
        if(length(mledist.par$data)==0){
            x.fit[[mm]] <- empty.fit ## if there are no data in the fitting period: NA
            x.fit.monitor[mm] <- 4
            if(warn){
                warning("all values in month ",mm," are 'NA'")
            }
        } else if(all(mledist.par$data==mledist.par$data[1])){ ## if all values in calibration period are eaqual...
            x.fit[[mm]] <- empty.fit
            x.fit.monitor[mm] <- 5
            if(warn){
                warning("all values in month ",mm," are constant, distribution not defined")
            }
        } else {
            ## a) fit distribution
            if(p0){
                np0 <- sum(mledist.par$data==0)
                nn <- length(mledist.par$data)
                if(!p0.center.mass){
                    p0.est <- np0/nn
                } else {
                    p0.est <- np0/(nn+1)                    
                }
                if(p0.est>0){
                    mledist.par$data <- mledist.par$data[mledist.par$data>0]
                    ## Catch that adds a single value close to zero if there are zeros in the
                    ## fitting period.  This forces the distribution to tend towards zero,
                    ## preventing a "gap" betwen 0 and data.
                    mledist.par$data <- c(mledist.par$data,0.01*min(mledist.par$data,na.rm=TRUE))
                }
            }
            mledist.par$start <- start.fun(x=mledist.par$data,distr=distr)
            fail.value <- mledist.par$start
            fail.value <- unlist(fail.value) 
            if(!start.fun.fix){
                fail.value[] <- NA
            }
            if(any(is.na(unlist(mledist.par$start)))){
                if(any(is.na(unlist(mledist.par$start))) & warn)
                    warning("starting values in month ",mm,
                            " not properly estimated \n parameter are NA.")
                fail.value[] <- NA
                x.fit[[mm]] <- fail.value
                x.fit.monitor[mm] <- 1
            } else {
                sink <- capture.output(
                    d.fit <- try(suppressWarnings(do.call(mledist,mledist.par)),silent=TRUE)
                    )
                ## b) save values
                if(class(d.fit)=="try-error"){
                    x.fit[[mm]] <- fail.value
                    x.fit.monitor[mm] <- 2
                } else if(d.fit$convergence > 0){
                    x.fit[[mm]] <- fail.value
                    x.fit.monitor[mm] <- 3
                } else {
                    x.fit[[mm]] <- d.fit$estimate
                }
                if(x.fit.monitor[mm] > 1 & start.fun.fix & warn)
                    warning("maximum likelihood estimation failed for month ", mm ,
                            "\n use starting values instead.")
                if(x.fit.monitor[mm] > 1 & !start.fun.fix & warn)
                    warning("maximum likelihood estimation failed for month ", mm ,
                            "\n parameters are set to NA.")
                if(p0){
                    if(!p0.center.mass){                   
                        x.fit[[mm]] <- c(x.fit[[mm]],P0=p0.est)
                    } else {
                        x.fit[[mm]] <- c(x.fit[[mm]],P0=p0.est,N.P0=np0,N=nn)
                    }
                }
                if(any(is.na(x.fit[[mm]]))){
                    x.fit[[mm]][] <- NA
                }
            }
        }
    }
    op <- list(dist.para=do.call(cbind,x.fit),
               dist.para.flag=x.fit.monitor,
               time.scale=time.scale,
               distr=distr,
               p0=p0,
               p0.center.mass=p0.center.mass,
               scaling=scale.val,
               call=match.call())
    class(op) <- "fitSCI"
    return(op)
}

transformSCI <- function(x,...)
    UseMethod("transformSCI")

transformSCI.default <- function(x,first.mon,obj,sci.limit=Inf,warn=TRUE,...){
    x <- as.numeric(x)
    x <- x/obj$scaling
    nn <- length(x)
    time.index <- 1:nn
    mmon <- .mon.vect(first.mon=first.mon,n=nn)
    ## 1. backward moving average
    time.scale <- obj$time.scale
    if(time.scale>1){
        ffilt <- rep(1,time.scale)/time.scale
        x <- as.numeric(filter(x,filter=ffilt,method="convolution",sides=1))
    }
    ## 2. estimate probability based on fitted distribution
    pdistr <- obj$distr
    pdistr <- match.fun(paste("p",pdistr,sep=""))
    if(obj$p0){
        npar <- nrow(obj$dist.para)
        PP0 <- obj$dist.para["P0",]
        if(!is.null(dim(PP0)))
            PP0 <- PP0[nrow(PP0),]
        if(!obj$p0.center.mass){
            distpar <- obj$dist.para[-npar,]
        } else {
            N.P0 <- obj$dist.para["N.P0",]
            if(!is.null(dim(N.P0)))
                N.P0 <- N.P0[nrow(N.P0),]
            NN <- obj$dist.para["N",]
            if(!is.null(dim(NN)))
                NN <- N.P0[nrow(NN),]
            distpar <- obj$dist.para[1:(npar-3),]
        }
    } else {
        distpar <- obj$dist.para
    }
    for(mm in 1:12){
        xm <- x[mmon==mm]
        if(any(is.na(distpar[,mm]))){
            xm[] <- NA
            if(warn)
                warning("Parameter for month ",mm," is NA")
        }  else {              
            xm <- try({
                xh <- do.call(pdistr,c(list(xm),distpar[,mm]))
                if(obj$p0){
                    xh <- PP0[mm] + (1-PP0[mm])*xh
                    if(obj$p0.center.mass){
                        ## browser()
                        xh[xm==0] <- (N.P0[mm]+1)/(2*(NN[mm]+1))                        
                    }
                }               
                xh
            },silent=TRUE)
            if(class(xm)=="try-error"){
                xm[] <- NA
                if(warn)
                    warning("transformation for month ", mm," faild")
            }
        }          
        x[mmon==mm] <- xm
    }
    ## 3. transform to "normal"
    x <- qnorm(x)
    ## 4. apply truncation
    x[x>sci.limit] <- sci.limit
    x[x< -sci.limit] <- -sci.limit
    return(x)
}
