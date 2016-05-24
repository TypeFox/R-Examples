## d <- twinsim(5000,b1=c(0),acde=c(0.6,0.2,0,0.2)); d$group <- factor(d$y0,labels=c("A","B"))
## b <- twinlm.strata(y~1,data=d,id="id",zyg="zyg",DZ="DZ",var="x1",quantiles=c(0,0.25,0.5,0.75,1))
## b2 <- twinlm.strata(y~1,data=d,id="id",zyg="zyg",DZ="DZ",var="group",type="ae")
## plot(b2,which=c(2,3),ylim=c(0,1),col=c("darkblue","darkred"),lwd=3,delta=.02,pch=16)
## b2 <- twinlm.strata(y0~1,data=d,id="id",zyg="zyg",DZ="DZ",var="x1",quantiles=c(0,0.25,0.5,0.75),binary=TRUE)
## plot(b,which=5:6,ylim=c(0,1),type="l")
## plot(b,which=5:6,ylim=c(0,1),col=c("darkred","darkblue"),legend=c("MZ","DZ"),lty=1:2)

##' @export
coef.timemets <- function(object,...) {    
    res <- unlist(lapply(object$summary,function(x) x$estimate[,1]))
    if (!is.null(names(object$var))) {
        nn <- names(object$var)[seq(length(object$summary))]
    } else {
        nn <- object$var[seq(length(object$summary))]
    }
    names(res) <- paste(names(res),rep(nn,length(res)/length(nn)),sep=":")
    return(res)
}

##' @export
vcov.timemets <- function(object,...) {
    if (object$type=="strata") {
        return(do.call(blockdiag,
                       lapply(object$summary, function(x) x$vcov)))
    }
}


##' @export
twinlm.strata <- function(formula,data,var,breaks,quantiles,...) {
    varname <- ""
    if (is.character(var) && length(var)==1) {
        varname <- var
        var <- data[,var]
    }
    ## 
    if (!missing(quantiles)) breaks <- quantile(var,quantiles,include.lowest=TRUE)
    if (!missing(breaks)) {
        var <- cut(var,c(breaks))
    }
    if (!inherits(var,c("character","factor","logical"))) stop("Expected factor, character, or logical vector")
    lev <- na.omit(unique(var))
    res <- list()
    for (i in seq_along(lev)) {
        tau <- lev[i]
        if (length(lev)>1) message(tau)
        data0 <- data[var==tau,,drop=FALSE]
        suppressWarnings(b <- twinlm(formula,data=data0,...))
        res <- c(res,list(summary(b)))
    }
    if (length(lev)==1) return(b)
    if (!missing(breaks)) lev <- breaks
    coef <- c(lapply(res,function(x) x$all),list(res[[length(res)]]$all))
    res <- list(varname=varname,var=lev,coef=coef,summary=res,type="strata")
    class(res) <- "timemets"
    return(res)
}


## data(prt)
## bb <- twinlm.time(cancer~country,data=prt,id="id",zyg="zyg",DZ="DZ",cens.formula=Surv(time,status==0)~zyg,breaks=seq(70,90,by=4))
## plot(bb,which=c(7,11),ylim=c(0,28),legendpos="topright",col=c("darkred","darkblue"),lty=c(1,2),legend=c("MZ","DZ"),ylab="Relative recurrence risk ratio")
## plot(bb,which=c(7,11),ylim=c(0,28),legendpos="topright",col=c("darkred","darkblue"),lty=c(1,2),legend=c("MZ","DZ"),ylab="Relative recurrence risk ratio",type="l")

##' @export
twinlm.time <- function(formula,...) {
    biprobit.time(formula,estimator="bptwin",...)
}

##' @export
bptwin.time <- function(formula,...) {
    biprobit.time(formula,estimator="bptwin",...)
}

##' @export
summary.timemets <- function(object,which=seq(nrow(object$coef[[1]])),...) {
    res <- list()
    for (i in which) {    
        rr <- matrix(unlist(lapply(object$coef,function(z) z[i,])),ncol=3,byrow=TRUE)
        colnames(rr) <- colnames(object$coef[[1]])
        rr <- cbind(object$var,as.data.frame(rr[seq_along(object$var),,drop=FALSE]))
        colnames(rr)[1] <- object$varname
        res <- c(res,list(rr))
    }
    names(res) <- rownames(object$coef[[1]])[which]
    return(res)
}


##' @export
print.timemets <- function(x,tail,row.names=FALSE,digits=4,width=10,...) {
    res <- summary(x,...)
    if (!is.null(x$summary[[1]]$ncontrast) && x$summary[[1]]$ncontrast>1) {
        cat("Contrasts:\n")
        for (i in seq(x$summary[[1]]$ncontrasts)) {
            cat("   c",i,":\n",sep="")
            cat("\tDependence ", x$summary[[1]]$par[[i]]$corref, "\n")
            if (x$summary[[1]]$model$eqmarg) {
                cat("\tMean       ", x$summary[[1]]$par[[i]]$mref1, "\n")
            } else {
                cat("\tMean 1     ", x$summary[[1]]$par[[i]]$mref1, "\n")
                cat("\tMean 2     ", x$summary[[1]]$par[[i]]$mref2, "\n")
            }
        }
    }
    ## }
    ## if (mcontr2 || (x$contrast & !x$model$eqmarg)) {
    ##     cat("\tMean 1       ", x$par[[i]]$mref1, "\n")
    ##     cat("\tMean 2       ", x$par[[i]]$mref2, "\n")
    ## } 
    ## if (mcontr1 || (x$contrast & x$model$eqmarg))
    ##     cat("\tMean         ", x$par[[i]]$mref1, "\n")
    
    
    M <- res[[1]][,1]
    nn <- c()
    for (i in seq_along(res)) {
        nn <- c(nn,names(res)[i])
        M <- cbind(M,res[[i]][,2])
    }
    nn0 <- cbind(paste(seq(ncol(M)-1),":",nn,sep=""))
    colnames(nn0) <- ""; rownames(nn0) <- rep("",nrow(nn0))
    print(nn0,quote=FALSE)
    cat("\n")
    nn <- unlist(lapply(nn,
                        function(x) {
                            res <- toString(x,width)
                            if (nchar(res)==nchar(x)) return(x)
                            substr(res,1,nchar(res)-1)
                        }))
    nn <- paste(seq(ncol(M)-1),":",nn,sep="")
    colnames(M) <- c("Time",nn)
    if (!missing(tail)) {
        print(utils::tail(round(M,digits=digits),tail),row.names=row.names)
    } else {
        print(round(M,digits=digits),row.names=row.names)
    }
    invisible(M)
}

##' @export
plot.timemets <- function(x,...,which=1,
                            type="s",
                            lwd=2,lty=1,col,fillcol,alpha=0.2,
                            xlab=x$varname,
                            ylab="",idx=seq_along(x$var),
                            lasttick=TRUE,add=FALSE,
                            legend=TRUE,legendpos="topleft") {
    ss <- summary(x,which)
    if (missing(col)) col <- seq_along(which)
    if (length(col)==1) col <- rep(col,length(which))
    if (length(lwd)==1) lwd <- rep(lwd,length(which))
    if (length(lty)==1) lty <- rep(lty,length(which))
    if (alpha>0 & missing(fillcol)) fillcol <- Col(col,alpha)
    count <- 0
    if (add) dev <- devcoords()
    for (tt in seq_along(which)) {
        count <- count+1
        zz <- ss[[tt]][idx,,drop=FALSE]
        if (!add) {
            plot(zz[,1:2,drop=FALSE],lty=0,
                 ylab=ylab,xlab=xlab,type="n",...)
            dev <- devcoords()
        }
        if (lasttick && type=="s" && !is.factor(zz[,1])) {
            zz2 <- rbind(zz,tail(zz,1))
            zz2[nrow(zz2),1] <- dev$fig.x2
            confband(zz2[,1],zz2[,3],zz2[,4],polygon=TRUE,step=(type=="s"),col=fillcol[count],border=0)
        } else {
            if (is.factor(zz[,1])) {
                confband(x=seq(nrow(zz)),lower=zz[,3],upper=zz[,4],center=zz[,2],col=col[count],lwd=lwd[count],lty=lty[count],...)
            } else {
                confband(zz[,1],zz[,3],zz[,4],polygon=TRUE,step=(type=="s"),col=fillcol[count],border=0)
            }
        }
        add <- TRUE
        ## if (lasttick && type=="s") {
        ##     axis(4,at=zz[nrow(zz),3],labels=FALSE,
        ##          lwd=lwd[count],col=fillcol[count],tcl=.5,...)
        ##     axis(4,at=zz[nrow(zz),4],labels=FALSE,
        ##          lwd=lwd[count],col=fillcol[count],tcl=.5,...)
        ## }
        if (!is.factor(zz[,1]))
            lines(zz[,1:2,drop=FALSE],lwd=lwd[count],lty=lty[count],col=col[count],type=type,...)        
    }
    if (!is.null(legend) || (is.logical(legend) && !legend[1])) {
        if (is.logical(legend) || length(legend)==1) legend <- rownames(x$coef[[1]])[which]
        graphics::legend(legendpos,legend=legend,col=col,lwd=lwd,lty=lty)
    }
    invisible(x)    
}


##' @export
bootstrap.timemets <- function(x,R=1000,...) {
    
}
