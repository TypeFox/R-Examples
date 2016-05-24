##' Wrapper function for mclapply
##'
##' @export
##' @param x function or 'sim' object
##' @param R Number of replications
##' @param f Optional function (i.e., if x is a matrix)
##' @param colnames Optional column names
##' @param messages Messages
##' @param mc.cores Number of cores to use
##' @param cl (optional) cluster to use for parallelization
##' @param blocksize Split computations in blocks
##' @param type type=0 is an alias for messages=1,mc.cores=1,blocksize=R
##' @param seed (optional) Seed (needed with cl=TRUE)
##' @param ... Additional arguments to (mc)mapply
##' @aliases sim.default summary.sim
##' @examples
##' m <- lvm(y~x+e)
##' distribution(m,~y) <- 0
##' distribution(m,~x) <- uniform.lvm(a=-1.1,b=1.1)
##' transform(m,e~x) <- function(x) (1*x^4)*rnorm(length(x),sd=1)
##'
##' onerun <- function(iter=NULL,...,n=2e3,b0=1,idx=2) {
##'     d <- sim(m,n,p=c("y~x"=b0))
##'     l <- lm(y~x,d)
##'     res <- c(coef(summary(l))[idx,1:2],
##'              confint(l)[idx,],
##'              estimate(l,only.coef=TRUE)[idx,2:4])
##'     names(res) <- c("Estimate","Model.se","Model.lo","Model.hi",
##'                     "Sandwich.se","Sandwich.lo","Sandwich.hi")
##'     res
##' }
##'
##' val <- sim(onerun,R=10,b0=1,messages=0,mc.cores=1)
##' val
##' val <- sim(val,R=40,b0=1,mc.cores=1) ## append results
##'
##' summary(val,estimate=c(1,1),confint=c(3,4,6,7),true=c(1,1))
##' summary(val,estimate=c(1,1),se=c(2,5),names=c("Model","Sandwich"))
##'
##' if (interactive()) {
##'     plot(val,estimate=1,c(2,5),true=1,names=c("Model","Sandwich"),polygon=FALSE)
##'     plot(val,estimate=c(1,1),se=c(2,5),main=NULL,
##'          true=c(1,1),names=c("Model","Sandwich"),
##'          line.lwd=1,density.col=c("gray20","gray60"),
##'          rug=FALSE)
##'     plot(val,estimate=c(1,1),se=c(2,5),true=c(1,1),
##'          names=c("Model","Sandwich"))
##' }
sim.default <- function(x=NULL,R=100,f=NULL,colnames=NULL,messages=1L,mc.cores,blocksize=2L*mc.cores,cl,type=1L,seed=NULL,...) {    
    if (missing(mc.cores) || .Platform$OS.type=="windows") {
        if (.Platform$OS.type=="windows") { ## Disable parallel processing on windows
            mc.cores <- 1L
        } else {
            mc.cores <- getOption("mc.cores",parallel::detectCores())
        }
    }
    if (type==0L) {
        mc.cores <- 1L
        if (inherits(R,c("matrix","data.frame")) || length(R)>1) {
            blocksize <- NROW(R)
        } else {
            blocksize <- R
        }
        messages <- 0
    }
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        runif(1)
    if (is.null(seed))
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    if (mc.cores>1L || !missing(cl)) requireNamespace("parallel",quietly=TRUE)
    newcl <- FALSE
    if (!missing(cl) && is.logical(cl)) {
        if (.Platform$OS.type=="windows" || TRUE) { ## Don't fork processes on windows
            cl <- NULL
            mc.cores <- 1
        } else {
            if (cl) {
                cl <- parallel::makeForkCluster(mc.cores)
                if (!is.null(seed)) parallel::clusterSetRNGStream(cl,seed)
                newcl <- TRUE
            }
        }
    }
    olddata <- NULL
    dots <- list(...)
    mycall <- match.call(expand.dots=FALSE)
    if (inherits(x,c("data.frame","matrix"))) olddata <- x
    if (inherits(x,"sim")) {        
        oldcall <- attr(x,"call")
        x <- attr(x,"f")
        if (!is.null(f)) x <- f
        ex <- oldcall[["..."]]
        for (nn in setdiff(names(ex),names(dots))) {
            dots[[nn]] <- ex[[nn]]
            val <- list(ex[[nn]]); names(val) <- nn
            mycall[["..."]] <- c(mycall[["..."]],list(val))
        }
        
    } else {
        if (!is.null(f)) x <- f
        if (!is.function(x)) stop("Expected a function or 'sim' object.")
    }
    if (is.null(x)) stop("Must give new function argument 'f'.")
    res <- val <- NULL
    on.exit({
        if (messages>0) close(pb)
        if (newcl) parallel::stopCluster(cl)
        if (is.null(colnames) && !is.null(val)) {
            if (is.matrix(val[[1]])) {
                colnames <- base::colnames(val[[1]])
            } else {
                colnames <- names(val[[1]])
            }
        }
        base::colnames(res) <- colnames
        if (!is.null(olddata)) res <- rbind(olddata,res)
        attr(res,"call") <- mycall
        attr(res,"f") <- x
        class(res) <- c("sim","matrix")
        if (idx.done<R) {
            res <- res[seq(idx.done),,drop=FALSE]
        }
        return(res)
    })
    if (inherits(R,c("matrix","data.frame")) || length(R)>1) {
        parval <- as.data.frame(R)
        if (is.vector(R)) names(parval) <- NULL
        else if (inherits(R,c("matrix","data.frame"))) names(parval) <- colnames(R)
        R <- NROW(parval)        
    } else {
        parval <- as.data.frame(1:R)
        names(parval) <- NULL
    }
    nfolds <- max(1,round(R/blocksize))
    idx <- split(1:R,sort((1:R)%%nfolds))
    idx.done <- 0
    count <- 0
    if (messages>0) pb <- txtProgressBar(style=lava.options()$progressbarstyle,width=40)
    for (ii in idx) {
        count <- count+1
        if (!missing(cl) && !is.null(cl)) {
            pp <- c(as.list(parval[ii,,drop=FALSE]),dots,list(cl=cl,fun=x,SIMPLIFY=FALSE))
        } else {
            pp <- c(as.list(parval[ii,,drop=FALSE]),dots,list(mc.cores=mc.cores,FUN=x,SIMPLIFY=FALSE))
        }
        if (mc.cores>1) {
            if (!missing(cl) && !is.null(cl)) {
                val <- do.call(parallel::clusterMap,pp)
            } else {
                val <- do.call(parallel::mcmapply,pp)
            }
        } else {
            val <- do.call(mapply,pp)
        }
        if (messages>0) 
            setTxtProgressBar(pb, count/length(idx))
        if (is.null(res)) {
            ##res <- array(NA,dim=c(R,dim(val[[1]])),dimnames=c(list(NULL),dimnames(val[[1]]),NULL))
            res <- matrix(NA,ncol=length(val[[1]]),nrow=R)
        }
        res[ii,] <- Reduce(rbind,val)
        ##rr <- abind::abind(val,along=length(dim(res)))
        ##res[ii,] <- abind(val,along=length(dim(res)))
        idx.done <- max(ii)
    }
}

##' @export
"[.sim" <- function (x, i, j, drop = FALSE) {
    atr <- attributes(x)
    class(x) <- "matrix"
    x <- NextMethod("[",drop=FALSE)
    atr.keep <- "call"
    if (missing(j)) atr.keep <- c(atr.keep,"f")
    attributes(x)[atr.keep] <- atr[atr.keep]
    class(x) <- c("sim","matrix")
    x
}

##' @export
print.sim <- function(x,...) {
    attr(x,"f") <- attr(x,"call") <- NULL
    class(x) <- "matrix"
    print(x,...)
}



##' Plot sim object
##'
##' @examples
##' n <- 1000
##' val <- cbind(est1=rnorm(n,sd=1),est2=rnorm(n,sd=0.2),est3=rnorm(n,1,sd=0.5),
##'              sd1=runif(n,0.8,1.2),sd2=runif(n,0.1,0.3),sd3=runif(n,0.25,0.75))
##'
##' plot.sim(val,estimate=c(1,2),true=c(0,0),se=c(4,5),equal=TRUE)
##' plot.sim(val,estimate=c(1,3),true=c(0,1),se=c(4,6),density.xlim=c(-3,3),ylim=c(-3,3))
##' plot.sim(val,estimate=c(1,2),true=c(0,0),se=c(4,5),equal=TRUE,plot.type="single")
##' plot.sim(val,estimate=c(1),se=c(4,5,6),plot.type="single")
##' plot.sim(val,estimate=c(1,2,3),equal=TRUE)
##' plot.sim(val,estimate=c(1,2,3),equal=TRUE,byrow=TRUE)
##' plot.sim(val,estimate=c(1,2,3),plot.type="single")
##' plot.sim(val,estimate=1,se=c(3,4,5),plot.type="single")
##'
##' density.sim(val,estimate=c(1,2,3))
##' @param x sim object
##' @param plot.type Single or multiple plots
##' @param ... Additional graphical arguments
##' @aliases density.sim plot.sim
##' @export
##' @export density.sim
density.sim <- function(x,plot.type="single",...) {
    plot.sim(x,...,scatter.plot=FALSE,plot.type=plot.type)
}

##' @export
##' @export plot.sim
plot.sim <- function(x,estimate,se=NULL,true=NULL,
                     names=NULL,
                     auto.layout=TRUE,
                     byrow=FALSE,
                     type="p",
                     ask=grDevices::dev.interactive(),
                     line.col=1, line.lwd=1.8,
                     col=c("gray60","orange","darkblue","seagreen","darkred"),
                     pch=16,cex=0.5,lty=1,
                     true.lty=2,true.col="gray70",true.lwd=1.2,
                     legend,
                     legendpos="topleft",
                     cex.legend=0.6,
                     plot.type=c("multiple","single"),
                     polygon=TRUE,
                     cex.axis=0.5,
                     alpha=0.5,
                     rug=TRUE,
                     rug.alpha=0.5,                     
                     main,
                     cex.main=1,
                     equal=FALSE,
                     delta=1.15,
                     ylim=NULL,
                     ylab="Estimate",
                     density.ylab="Density",
                     density.ylim=NULL,
                     density.xlim=NULL,
                     density.plot=TRUE,
                     scatter.plot=TRUE,
                     running.mean=scatter.plot,
                     density.alpha=0.2,
                     border=density.col,
                     density.lty,
                     density.col=col,
                     density.lwd=0.4,
                     xlab="",...) {

    if (missing(estimate)) {
        estimate <- seq(ncol(x))
    }
    if (is.null(estimate)) {
        av <- apply(x[,drop=FALSE],2,function(z) cumsum(z)/seq(length(z)))
        graphics::matplot(x,type="p",pch=pch, cex=cex, col=col,...)
        graphics::matlines(av,type="l",col=col,lty=lty,...)
        if (!is.null(true)) abline(h=true,lty=true.lty,...)
        if (missing(legend)) legend <- colnames(x)
        if (!is.null(legend))
            graphics::legend(legendpos,legend=legend,bg="white",col=col,lty=lty,pch=pch,...)
        return(invisible(NULL))
    }
    if (is.character(estimate)) {
        estimate <- match(estimate,colnames(x))
    }    

    K <- length(estimate)
    est <- tru <- c()
    if (length(se)>0) {
        if (K==1 && !is.list(se))
            se <- list(se)
        else se <- as.list(se)
    } else {
        est <- estimate; tru <- true
    }
    for (i in seq_along(estimate)) {
        est <- c(est,list(rep(estimate[i],length(se[[i]]))))
        if (!is.null(true)) tru <- c(tru,list(rep(true[i],length(se[[i]]))))
    }
    
    if (length(se)>0) {
        for (i in seq_along(se)) {
            if (is.character(se[[i]])) se[[i]] <- match(se[[i]],colnames(x))
        }
    }
    
    ss <- summary.sim(x,estimate=unlist(est),se=unlist(se),true=unlist(tru),names=names)
    oldpar <- NULL
    on.exit({
        par(oldpar)
        return(invisible(ss))
    })

    single <- tolower(plot.type[1])=="single"

    if (auto.layout) {
        nc <- (scatter.plot || running.mean) + density.plot
        nr <- min(6,K)
        if (single) nr <- 1
        oma.multi = c(2, 0, 2, 0)
        mar.multi = c(1.5, 4.1, 1, 1)
        oldpar <- par(mar=mar.multi, oma=oma.multi,
                      cex.axis=cex.axis,las=1,
                      ask=FALSE)
        if (byrow) {
            par(mfrow=c(nr,nc))
        } else {
            par(mfcol=c(nc,nr))
        }
    }

    dys <- c()
    maxdy <- 0
    if (density.plot)
        for (i in seq(K)) {
            ii <- estimate[i]
            y <- as.vector(x[,ii])
            dy <- stats::density(y)
            dys <- c(dys,list(dy))
            maxdy <- max(maxdy,dy$y)
        }

    if (equal || single) {
        if (is.null(ylim)) {
            rg <- range(x[,estimate])
            rg <- rg+c(-1,1)*abs(diff(rg)*(delta-1))
            ylim <-  rep(list(rg),K)
        }
        if (density.plot) {
            if (is.null(density.xlim)) density.xlim <- ylim
            if (is.null(density.ylim)) density.ylim <- rep(list(c(0,maxdy*delta)),K)
        }
    }

    if (!is.null(ylim)) {
        if (!is.list(ylim)) ylim <- list(ylim)
        ylim <- rep(ylim,length.out=K)
    }
    ylab <- rep(ylab,length.out=K)
    if (!is.null(density.ylim)) {
        if (!is.list(density.ylim)) density.ylim <- list(density.ylim)
        density.ylim <- rep(density.ylim,length.out=K)
    }
    if (!is.null(density.xlim)) {
        if (!is.list(density.xlim)) density.xlim <- list(density.xlim)
        density.xlim <- rep(density.xlim,length.out=K)
    }
    if (missing(main)) {
        main <- NULL
        if (K>1 && !single) main <- colnames(ss)
    }
    if (!is.null(main)) main <- rep(main,length.out=K)
    if (missing(density.lty)) {
        density.lty <- rep(1,K)
        if (single || !polygon) {
            density.lty <- 1:20
        }
    }

    my.scatter.sim <- function(i,add=FALSE,colors,...) {
        ii <- estimate[i]
        if (!missing(colors)) {
            col <- line.col <- true.col <- colors[1]
        }
        y <- as.vector(x[,ii])
        args <- list(y,ylab=ylab[i],col=Col(col[1],alpha),cex=cex,pch=pch,type=type)
        if (!is.null(ylim)) args <- c(args,list(ylim=ylim[[i]]))
        if (scatter.plot) {
            if (!add) {
                do.call(graphics::plot,args)
            } else {
                do.call(graphics::points,args)
            }
        }
        if (running.mean) {
            lines(cumsum(y)/seq_along(y),col=line.col[1],lwd=line.lwd,lty=lty)
            if (!is.null(true))
                abline(h=true[i],lty=true.lty,col=true.col[1],lwd=true.lwd)
        }
    }

    my.density.sim <- function(i,add=FALSE,colors,alphas=density.alpha,auto.legend=TRUE,...) {
        ii <- estimate[i]
        y <- as.vector(x[,ii])
        if (!missing(colors)) {
            density.col <- border <- colors
            col <- true.col <- colors
        }
        if (density.plot) {
            dy <- stats::density(y)
            if (is.null(density.ylim)) {
                density.ylim0 <- c(0,max(dy$y)*delta)
            } else {
                density.ylim0 <- density.ylim[[i]]
            }
            if (is.null(density.xlim)) {
                density.xlim0 <- range(dy$x)
            } else {
                density.xlim0 <- density.xlim[[i]]
            }
            if (!add) graphics::plot(0,0,type="n",main="",ylab=density.ylab,xlab=xlab,ylim=density.ylim0,xlim=density.xlim0)
            if (polygon) {
                with(dy, graphics::polygon(c(x,rev(x)),c(y,rep(0,length(y))),col=Col(density.col[1],alpha=alphas[1]),border=NA))
                if (!is.null(border)) with(dy, lines(x,y,col=border[1],lty=density.lty[1],lwd=density.lwd[1]))
            } else {
                graphics::lines(dy,main="",lty=density.lty[1],col=density.col[1],lwd=density.lwd[1])
            }
            if (rug) graphics::rug(y,col=Col(col[1],rug.alpha[1]))
            if (!is.null(main) && !(running.mean || scatter.plot)) {
                title(main[i],cex.main=cex.main)
            }
            if (!is.null(true)) {
                abline(v=true[i],lty=true.lty,col=true.col,lwd=true.lwd)
            }

            if (!is.null(se)) {
                se.pos <- match(se[[i]],unlist(se))
                ns <- length(se.pos)+1
                se.alpha <- rep(alphas,length.out=ns)[-1]
                se.border <- rep(border,length.out=ns)[-1]
                se.col <- rep(density.col,length.out=ns)[-1]
                se.lty <- rep(density.lty,length.out=ns)[-1]
                se.lwd <- rep(density.lwd,length.out=ns)[-1]
                xx <- dy$x
                 for (j in seq_along(se.pos)) {
                    if (polygon) {
                        yy <- dnorm(xx,mean=ss["Mean",i],sd=ss["SE",se.pos[j]])
                        if (se.alpha[j]>0) graphics::polygon(c(xx,rev(xx)),c(yy,rep(0,length(yy))),col=Col(se.col[j],alpha=se.alpha[j]),border=NA)
                        if (!is.null(border)) lines(xx,yy,col=se.border[j],lty=se.lty[j],lwd=se.lwd[j])
                    } else {
                        graphics::curve(dnorm(x,mean=ss["Mean",i],sd=ss["SE",se.pos[j]]),lwd=se.lwd[j],lty=se.lty[j],col=se.col[j],add=TRUE)
                    }
                 }
                if (auto.legend) legend <- c("Kernel",colnames(ss)[se.pos])
                if (!is.null(legend)) {
                    if (polygon) {
                        dcol <- c(density.col[1],se.col)
                        graphics::legend(legendpos,legend,
                                         fill=Col(dcol,density.alpha),border=dcol,cex=cex.legend)
                    } else {
                        graphics::legend(legendpos,legend,
                                         col=c(density.col[1],se.col),
                                         lty=c(density.lty[1],se.lty),
                                         lwd=c(density.lwd[1],se.lwd),cex=cex.legend)
                    }
                }
            }

        }
    }

    if (single) {
        N <- K
        nk <- lapply(se,length)
        if (!is.null(se)) N <- sum(unlist(nk))+K
        col <- rep(col,length.out=K)
        for (i in seq(K)) {
            my.scatter.sim(i,add=(i>1),colors=col[i])
        }
        if (!is.null(main) && !byrow) {
            title(main[1],cex.main=cex.main)
        }
        if (missing(legend)) legend <- colnames(x)[estimate]
        legendold <- legend
        legend <- NULL
        density.alpha <- rep(density.alpha,length.out=K)
        for (i in seq(K)) {
            alphas <- density.alpha[i]
            if (length(se)>0) alphas <- c(alphas,rep(0,nk[i]))
            my.density.sim(i,add=(i>1),colors=col[i],alphas=alphas,auto.legend=FALSE)
        }
        if (!is.null(legendold)) {
            legend <- rep(legendold,length.out=K)
            graphics::legend(legendpos,legend,
                             fill=Col(col,density.alpha),border=col,cex=cex.legend)
        }

    } else {
        for (i in seq(K)) {
            my.scatter.sim(i)
            if (!is.null(main) && !byrow && scatter.plot) {
                title(main[i],cex.main=cex.main)
            }
            my.density.sim(i,auto.legend=missing(legend))
            if (i==1 && ask) par(ask=ask)
        }
    }

}

##' @export
##' @export summary.sim
summary.sim <- function(object,estimate=NULL,se=NULL,confint=NULL,true=NULL,fun,names=NULL,...) {
    if (missing(fun)) fun <- function(x) {
        pp <- c(.025,.5,.975)
        res <- c(mean(x,na.rm=TRUE),sd(x,na.rm=TRUE),quantile(x,c(0,pp,1),na.rm=TRUE),
                 mean(is.na(x)))
        names(res) <- c("Mean","SD","Min",paste0(pp*100,"%"),"Max","Missing")
        res
    }
    if (is.null(estimate) && is.null(confint)) return(apply(object,2,fun))

    if (is.character(estimate)) {
        estimate <- match(estimate,colnames(object))
    }
    
    est <- apply(object[,estimate,drop=FALSE],2,
                 function(x) c(Mean=mean(x,na.rm=TRUE),Missing=mean(is.na(x)),SD=sd(x,na.rm=TRUE)))
    if (!is.null(true)) {
        if (length(true)!=length(estimate)) stop("'true' should be of same length as 'estimate'.")
        est <- rbind(rbind(True=true),rbind(Bias=true-est["Mean",]),
                     rbind(RMSE=((true-est["Mean",])^2+(est["SD",])^2)^.5),
                     est)
    }
    if (!is.null(se)) {
        if (is.character(se)) {
            se <- match(se,colnames(object))
        }
        if (length(se)!=length(estimate)) stop("'se' should be of same length as 'estimate'.")
        est <- rbind(est, SE=apply(object[,se,drop=FALSE],2,
                                  function(x) c(mean(x,na.rm=TRUE))))
        est <- rbind(est,"SE/SD"=est["SE",]/est["SD",])

    }
    if (!is.null(confint)) {
        if (is.character(confint)) {
            confint <- match(confint,colnames(object))
        }
        if (length(confint)!=2*length(estimate)) stop("'confint' should be of length 2*length(estimate).")
        Coverage <- c()
        for (i in seq_along(estimate)) {
            Coverage <- c(Coverage,
                          mean((object[,confint[2*(i-1)+1]]<true[i]) & (object[,confint[2*i]]>true[i]),na.rm=TRUE))
        }
        est <- rbind(est,Coverage=Coverage)
    }
    if (!is.null(names)) colnames(est) <- names
    colnames(est) <- make.unique(colnames(est))
    return(est)
}
