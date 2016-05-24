##' @export
`lifetable` <- function(x,...) UseMethod("lifetable")

##' Create simple life table 
##' 
##' @title Life table
##' @param x time formula (Surv) or matrix/data.frame with columns time,status or entry,exit,status
##' @param strata Strata
##' @param data data.frame
##' @param breaks Time intervals
##' @param confint If TRUE 95\% confidence limits are calculated
##' @param ... Additional arguments to lower level functions
##' @author Klaus K. Holst
##' @aliases lifetable lifetable.matrix lifetable.formula
##' @usage
##'  \method{lifetable}{matrix}(x, strata = list(), breaks = c(),
##'    confint = FALSE, ...)
##'
##'  \method{lifetable}{formula}(x, data=parent.frame(), breaks = c(),
##'    confint = FALSE, ...)
##' @examples
##' library(timereg)
##' data(TRACE)
##' 
##' d <- with(TRACE,lifetable(Surv(time,status==9)~sex+vf,breaks=c(0,0.2,0.5,8.5)))
##' summary(glm(events ~ offset(log(atrisk))+factor(int.end)*vf + sex*vf,
##'             data=d,poisson))
##' @export
lifetable.matrix <- function(x,strata=list(),breaks=c(),confint=FALSE,...) {
    if (ncol(x)==3) {
        status <- x[,3]
        entry <- x[,1]
        time <- x[,2]
    } else {
        status <- x[,2]
        time <- x[,1]
        entry <- rep(0,length(time))
    }
    LifeTable(time,status,entry,strata,breaks,confint,...)
}

##' @export
lifetable.formula <- function(x,data=parent.frame(),breaks=c(),confint=FALSE,...) {
    cl <- match.call()
    mf <- model.frame(x,data)
    Y <- model.extract(mf, "response")
    Terms <- terms(x, data = data)
    if (!is.Surv(Y)) stop("Expected a 'Surv'-object")
    if (ncol(Y)==2) {
        exit <- Y[,1]
        entry <- NULL ## rep(0,nrow(Y))
        status <- Y[,2]
    } else {
        entry <- Y[,1]
        exit <- Y[,2]
        status <- Y[,3]
    }
    strata <- list()
    X <- model.matrix(Terms, data)
    if (!is.null(intpos  <- attributes(Terms)$intercept))
        X <- X[,-intpos,drop=FALSE]
    if (ncol(X)>0) {
        strata <- as.list(model.frame(Terms,data)[,-1,drop=FALSE])
    }
    LifeTable(exit,status,entry,strata,breaks,confint,...)       
}

## lifetable.data.table <- function(x,entry,exit,status,strata,breaks,...) {
##     requireNamespace("data.table")
##     breaks <- sort(unique(breaks))
##     nbreaks <- length(breaks)
##     ## prs <- parse(text=paste0("pmin(x,",exit,")-pmax(x-1,",entry,")"))
##     ## dt[,eval(prs),by=strata]
##     system.time(dur1 <- x[,
##                            lapply(breaks[-1],function(x)
##                               eval(parse(text=paste0("pmin(x,",exit,")-pmax(x-1,",entry,")")))
##                                  ),by=strata])
##     dur1[dur1<0] <- NA
##     system.time(endur1 <- x[,
##                             lapply(breaks[-nbreaks],function(x)
##                                 eval(parse(text=paste0("x+1-pmax(x,",entry,")")))
##                                    ),by=strata])
##     enter <- dur1[,lapply(.SD,function(x) sum(!is.na(x))),by=strata]
##     atrisk <- dur1[,lapply(.SD,function(x) sum(x,na.rm=TRUE)),by=strata]
##     vidx <- seq(length(strata)+1,ncol(atrisk))
##     suppressWarnings(names(atrisk)[vidx] <- paste0("_R",seq_along(vidx)))   
##     system.time(eventcens <- x[,
##                                lapply(breaks[-1],function(x)
##                                    eval(parse(text=paste0("((pmin(x,",exit,")-pmax(x-1,",entry,"))<(x-pmax(x-1,",entry,")))*(1+",status,")")))
##                                            ),by=strata])
##     lost <- eventcens[,lapply(.SD,function(x) sum(x==1,na.rm=TRUE)),by=strata]
##     events <- eventcens[,lapply(.SD,function(x) sum(x==2,na.rm=TRUE)),by=strata]
##     suppressWarnings(names(events)[vidx] <- paste0("_E",seq_along(vidx)))
##     res <- fast.reshape(cbind(events,atrisk[,vidx,with=FALSE]))
##     return(res)
##     res <- subset(data.frame(enter=enter,
##                              atrisk=atrisk,
##                              lost=lost,
##                              events=events,
##                              ## int.start=c(-Inf,breaks),
##                              ## int.end=c(breaks,Inf),
##                              int.start=breaks[-length(breaks)],
##                              int.end=breaks[-1],
##                              surv=0,
##                              rate=events/atrisk))
## }


LifeTable <- function(time,status,entry=NULL,strata=list(),breaks=c(),confint=FALSE,interval=TRUE,mesg=FALSE) {    
    if (is.null(entry)) entry <- rep(0,NROW(time))
    if (mesg) message(dim(time))
    if ((is.matrix(time) || is.data.frame(time)) && ncol(time)>1) {
        if (ncol(time)==3) {
            status <- time[,3]
            entry <- time[,1]
            time <- time[,2]
        } else {
            status <- time[,2]
            time <- time[,1]
            entry <- rep(0,length(time))
        }
    }
    if (length(strata)>0) {
        a <- by(cbind(entry,time,status), strata,
                FUN=LifeTable, breaks=breaks, confint=confint)
        cl <- lapply(strata,class)
        nulls <- which(unlist(lapply(a,is.null)))
        nonnulls <- setdiff(seq_along(a),nulls)
        nn <- do.call("expand.grid",attributes(a)$dimnames)
        if (length(nulls)>0) nn <- nn[-nulls,,drop=FALSE]
        nam <- nn[rep(seq(NROW(nn)),each=NROW(a[[nonnulls[1]]])),,drop=FALSE]
        xx <- list()
        for (i in seq(ncol(nam))) {
            if (cl[i]%in%c("numeric","integer"))
                xx <- c(xx,list(as.numeric(as.character(nam[,i]))))
            else                    
                xx <- c(xx, list(do.call(paste("as.",as.character(cl[i]),sep=""),list(nam[,i]))))
        }
        xx <- as.data.frame(xx); colnames(xx) <- colnames(nam)
        res <- Reduce("rbind",a)
        res <- cbind(res,xx)
        return(res)
    }
    if (length(breaks)==0) breaks <- c(0,max(time,na.rm=TRUE))
    if (length(breaks)==1) breaks <- c(0,breaks)
    breaks <- sort(unique(breaks))
    ## en <- matrix(unlist(lapply(c(-Inf, breaks),function(x) pmax(x,entry))),
    ##               ncol=length(breaks)+1)
    ## ex <- matrix(unlist(lapply(c(breaks, Inf),function(x) pmin(x,time))),
    ##               ncol=length(breaks)+1)
    en <- matrix(unlist(lapply(breaks[-length(breaks)],function(x) pmax(x,entry))),
                  ncol=length(breaks)-1)
    ex <- matrix(unlist(lapply(breaks[-1],function(x) pmin(x,time))),
                  ncol=length(breaks)-1)
    dur <- ex-en
    endur <- rbind(breaks[-1])%x%cbind(rep(1,nrow(en)))-en
    ##endur <- rbind(c(breaks,Inf))%x%cbind(rep(1,nrow(en)))-en
    dur[dur<=0] <- NA
    enter <- colSums(!is.na(dur))
    atrisk <- colSums(dur,na.rm=TRUE)
    ##eventcens <- dur<endur
    eventcens <- rbind(apply(dur<endur,2,function(x) x*(status+1)))
    lost <- colSums(eventcens==1,na.rm=TRUE)
    events <- colSums(eventcens==2,na.rm=TRUE)
    rate <- events/atrisk; rate[is.nan(rate)] <- 0
    res <- subset(data.frame(enter=enter,
                             atrisk=atrisk,
                             lost=lost,
                             events=events,
                             ## int.start=c(-Inf,breaks),
                             ## int.end=c(breaks,Inf),
                             int.start=breaks[-length(breaks)],
                             int.end=breaks[-1],
                             rate=rate))
    if (interval) res$interval <- factor(paste0("[",res$int.start,";",res$int.end,")"))
    cumsum.na <- function(x,...) { x[is.na(x)] <- 0; cumsum(x) }
    if (length(strata)==0) res$surv <- with(res, exp(-cumsum.na(rate*(int.end-int.start))))    
    if (confint) {
        ff <- events ~ offset(log(atrisk))
        if (length(breaks)>2) ff <- update(ff,.~.+factor(int.end)-1)
        g <- glm(ff,data=res,poisson)
        suppressMessages(ci <- rbind(exp(stats::confint(g))))
        res[,"2.5%"] <- ci[,1]
        res[,"97.5%"] <- ci[,2]
    }
    res
}

eh <- function(formula,intervals,family=poisson(log),...) {
    if (missing(intervals)) stop("Please supply list of time-intervals")
    g <- glm(formula,family=family,...)
    structure(g,class=c("glm","eh"),intervals=intervals)
    
}

##' Summary for survival analyses via the 'lifetable' function
##'
##' Summary for survival analyses via the 'lifetable' function 
##' @title Extract survival estimates from lifetable analysis
##' @param object glm object (poisson regression)
##' @param ... Contrast arguments
##' @param timevar Name of time variable
##' @param time Time points (optional)
##' @param int.len Time interval length (optional)
##' @param confint If TRUE confidence limits are supplied
##' @param level Level of confidence limits
##' @param individual Individual predictions
##' @param length.out Length of time vector
##' @export
##' @author Klaus K. Holst
eventpois <- function(object,...,timevar,time,int.len,confint=FALSE,level=0.95,individual=FALSE,length.out=25) {      
    if (missing(timevar)) {
        timevar <- names(attributes(object)$intervals)[1]
        times <- attributes(object)$intervals
        int.len <- diff(times)
    }
    nn <- names(coef(object))
    if (is.null(timevar)) {
        idx <- unlist(sapply(c("\\(Intercept\\)","bs\\(","ns\\("),function(x) grep(x,nn)))
        keep <- setdiff(seq_along(nn),idx)
        cc <- estimate(object,exp,keep=keep,labels=nn[keep])$coefmat
        colnames(cc)[1] <- "RR"
        cc[,5] <- coef(summary(object))[keep,4]
        return(cc)
    }
    timevar_re0 <- gsub("\\$|\\^","",glob2rx(timevar))
    timevar_re <- paste0(timevar_re0,"[0-9]+\\.*[0-9]*")
    idx <- regexpr(timevar_re,nn)
    dots <- list(...)
    if (all(idx<0)) {
        timevar_re0 <- gsub("\\$|\\^","",glob2rx(paste0("factor(",timevar,")")))
        timevar_re <- paste0(timevar_re0,"[0-9]+\\.*[0-9]*")
        idx <- regexpr(timevar_re,nn)
    }
    
    tvar <- unique(regmatches(nn,idx))
    if (missing(time)) {        
        time <- sort(unique(as.numeric(gsub(timevar_re0,"",tvar))))
        if (length(time)==0) {
            i0 <- seq(nrow(object$data))
            ## browser()
            ## if ("rate"%in%colnames(object$data)) {
            ##     ii0 <- which(object$data[,"rate"]>0)
            ##     i0 <- intersect(i0,ii0)
            ## }
            ## for (i in seq_along(dots)) {
            ##      i0 <- intersect(i0,which(object$data[,names(dots)[i]]==dots[[i]]))
            ##  }
            rg <- range(object$data[i0,timevar],na.rm=TRUE)        
            time <- seq(rg[1],rg[2],length.out=length.out)
        } else {
            ##time <- seq(1,rg[2])
        }
    }
    if (missing(int.len)) {
        int.len <- diff(c(0,time))
    } else if (int.len==1) int.len <- rep(int.len,nrow(res))
    tt <- terms(object)
    offsetvar <- NULL
    if (attr(tt,"offset")) {
        offsetvar <- all.vars(tt)[attr(tt,"offset")]
        dots[[offsetvar]] <- 1
    }
    responsevar <- getoutcome(tt)
    vv <- setdiff(all.vars(formula(tt)),c(offsetvar,responsevar))
    dots0 <- dots
    dots[[timevar]] <- time
    dotnames <- c()
    if (individual) {
        newdata <- as.data.frame(dots)
        for (v in vv) {
            if (v%ni%names(dots)) newdata[,v] <- object$data[1,v]
            else dotnames <- c(dotnames,v)
        }    
    } else {
        for (v in vv) {            
            if (v%ni%names(dots)) {
                dots[[v]] <- object$data[1,v]
            } else {
                dotnames <- c(dotnames,v)
            }
        }        
        args <- c(list(`_data`=model.frame(object)),dots)
        newdata <- do.call(Expand, args)
        newdata <- dsort(newdata,c(setdiff(dotnames,timevar),timevar))
    }

    Terms <- delete.response(tt)    
    m <- model.frame(Terms, newdata, xlev = object$xlevels)
    X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
    ## NA-robust matrix product
    beta0 <- coef(object); beta0[is.na(beta0)] <- -.Machine$double.xmax    
    V <- vcov(object)    
    V0 <- structure(matrix(0,length(beta0),length(beta0)),dimnames=list(names(beta0),names(beta0)))
    idx <- which(rownames(V0)%in%rownames(V))
    V0[idx,idx] <- V
    coefs <- X%*%beta0
    res <- cbind(time-int.len,int.len,exp(coefs))
    colnames(res) <- c("time","int.len","rate")    
    if (nrow(coefs)<=length(time)) {
        if (!confint) {
            res <- cbind(res,cumsum(res[,3]*int.len))
            res <- cbind(res,exp(-res[,4]))
        } else {
            S0 <- X%*%V0%*%t(X)
            system.time(e0 <- estimate(NULL,coef=coefs,vcov=S0,function(x) cumsum(exp(x)*int.len))) ## Cumulative hazard
            system.time(e1 <- estimate(e0, function(x) exp(-x),level=level)) ## Survival
            res <- cbind(res,e0$coefmat[,1],e1$coefmat[,c(1,3:4)])
        }
        colnames(res)[4:5] <- c("chaz","surv") 
    }
    if (!is.null(offsetvar)) newdata[,offsetvar] <- NULL

    times <- list()
    if (!individual) {
        aa <- unique(newdata[,setdiff(dotnames,c(timevar,offsetvar)),drop=FALSE])
        vv <- colnames(aa)
        for (i in seq(nrow(aa))) {
            idx <- seq(nrow(object$data))
            for (j in seq_along(vv)) {
                val <- aa[i,j]
                if (is.factor(val)) val <- as.character(val)
                idx <- intersect(idx,which(object$data[,vv[j]]==val))
            }
            times <- c(times, list(sort(unique(object$data[idx,timevar]))))
        }
    }

    structure(as.data.frame(cbind(res,newdata)),
              args=dots0,
              variables=aa,
              times=times,
              individual=individual,
              class=c("eventpois","data.frame")
              )
}

##' @export
plot.eventpois <- function(x,var,confint=TRUE,cont=TRUE,length.out=200,type="surv",
                           line.type="l",lty=c(1:8),col=1,lwd=1,add=FALSE,
                           use.time=TRUE,
                           legend,                           
                           xlab="Time",ylab="Survival probability",xlim,ylim,...) {
    type <- agrep(type,c("clogclog","hazard","survival"))
    if (!add) {
        if (missing(xlim)) xlim <- c(0,max(x[,1]+x[,2]))
        if (missing(ylim)) {
            ylim <- switch(type,
                           '1'=c(-30,max(log(sum(x[,2]*x[,3])))),
                           '2'=c(0,max(x[,3])),
                           c(0,1))
        }
        plot(0,type="n",xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,...)
    }
    ##args <- attr(x,"args")
    ##aa <- expand.grid(args)
    aa <- attr(x,"variables")
    vv <- colnames(aa)    
    if (attr(x,"individual")) aa <- aa[1,,drop=FALSE]
    if (missing(legend)) legend <- nrow(aa)>1
    if (length(lty)<nrow(aa)) lty <- rep(lty,nrow(aa))
    if (length(col)<nrow(aa)) col <- rep(col,nrow(aa))
    if (length(lwd)<nrow(aa)) lwd <- rep(lwd,nrow(aa))
    for (i in seq(nrow(aa))) {
        idx <- seq(nrow(x))
        
        if (!attr(x,"individual"))
        for (j in seq_along(vv)) {
            val <- aa[i,j]
            if (is.factor(val)) val <- as.character(val)
            idx <- intersect(idx,which(x[,vv[j]]==val))
        }
        idx <- intersect(idx,which(x[,3]>0))
        e0 <- x[idx,,drop=FALSE]
        dt <- tail(e0,1)        
        time0 <- c(e0[,1],dt[,1]+dt[,2])
        if (use.time) {
            ##browser()
            time0 <- attr(x,"time")[[i]]-1
            e0 <- e0[na.omit(match(time0,e0[,"time"])),,drop=FALSE]
            time0 <- time0[match(e0[,"time"],time0)]
            dt <- tail(e0,1)
            time0 <- c(time0,dt[,1]+dt[,2])
        }
        chaz0 <- cumsum(e0[,3]*e0[,2])
        tt <- seq(min(time0),max(time0),length.out=length.out)
        ff <- approxfun(time0,c(0,chaz0),method="linear")
        if (type==3) {
            lines(tt,exp(-ff(tt)),lty=lty[i],col=col[i],lwd=lwd[i],...)
        } else if (type==1) {
            lines(tt,log(ff(tt)),lty=lty[i],col=col[i],lwd=lwd[i],...)
        } else {
            lines(time0,e0[,3],lty=lty[i],col=col[i],lwd=lwd[i],...)
        }
    }
    if (legend) {
        lab <- apply(aa,1,function(x) paste(x,collapse=","))
        graphics::legend("bottomleft",lab,lty=lty,col=col,lwd=lwd,bg="white")
    }
    return(aa)    
    ## if (surv) lines(x[,c(1,5)],lty=lty[1],type=type,...)
    ## else lines(x[,c(1,3)],lty=lty[1],type=type,...)
    ## if (ncol(x)>4 && confint && surv) {
    ##     lines(x[,c(1,6)],lty=lty[2],type=type,...)
    ##     lines(x[,c(1,7)],lty=lty[2],type=type,...)
    ## }
}


