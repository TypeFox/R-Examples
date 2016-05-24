##' @export
biprobit.time <- function(formula,data,id,...,
                          breaks=NULL,n.times=20,pairs.only=TRUE,fix.cens.weights=FALSE,
                          cens.formula,cens.model="aalen",weights="w",messages=FALSE,
                          return.data=FALSE,theta.formula=~1,trunc.weights="w2",
                          estimator="biprobit", summary.function) {
## {{{ 
    m <- match.call(expand.dots = FALSE)
    m <- m[match(c("","data"),names(m),nomatch = 0)]
    Terms <- terms(cens.formula,data=data)
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    M <- eval(m,envir=parent.frame())

    censtime <- model.extract(M, "response")
    if (pairs.only) {
        ii <- sort(as.matrix(na.omit(fast.reshape(seq(nrow(data)),id=data[,id]))))
        data <- data[ii,,drop=FALSE]
        censtime <- censtime[ii]
    }

    ltimes <- 0
    if (ncol(censtime)==3) {## {{{ ## Calculate probability of not being truncated via Clayton-Oakes (ad hoc combining causes)
        status <- censtime[,3]
        noncens <- !status
        time <- censtime[,2]
        ltimes <- censtime[,1]        
        data$truncsurv <- Surv(ltimes,time,noncens)
        trunc.formula <- update(formula,truncsurv~.)       
        ud.trunc <- aalen(trunc.formula,data=data,robust=0,n.sim=0,residuals=0,silent=1,max.clust=NULL,
                          clusters=data[,id], ...)
        X <- model.matrix(trunc.formula,data)
        ##dependX0 <- model.matrix(theta.formula,data)
        dependX0 <- X       
        twostage.fit <- two.stage(ud.trunc,
                                 data=data,robust=0,detail=0,
                                 theta.des=dependX0)#,Nit=20,step=1.0,notaylor=1)        
        Xnam <- colnames(X)
        ww <- fast.reshape(cbind(X,".num"=seq(nrow(X)),".lefttime"=ltimes),varying=c(".num",".lefttime"),id=data[,id])
        dependX <- as.matrix(ww[,Xnam,drop=FALSE])
        nottruncpair.prob <- predict.two.stage(twostage.fit,X=ww[,Xnam],
                                               times=ww[,".lefttime1"],times2=ww[,".lefttime2"],
                                               theta.des=dependX)$St1t2
        data[,trunc.weights] <- 0
        data[ww[,".num1"],trunc.weights] <- nottruncpair.prob
        data[ww[,".num2"],trunc.weights] <- nottruncpair.prob
        
    } else {  
        status <- censtime[,2]
        time <- censtime[,1]
    } ## }}} 

    outcome <- as.character(terms(formula)[[2]])
    jj <- jumptimes(time,data[,outcome],data[,id],sample=n.times)
    lastjump <- tail(jj,1)
    if (is.null(breaks)) {
        breaks <- jj
    }
    if (any(breaks>lastjump)) {
        breaks <- unique(pmin(breaks,tail(jj,1)))        
        if (messages) message("Looking at event before time ",max(breaks))
    }    

    outcome0 <- paste(outcome,"_dummy")
    res <- list(); k <- 0
    breaks <- rev(breaks)
    for (tau in breaks) {
        if (length(breaks)>1 && messages) message(tau)
        ## construct min(T_i,tau) or T_i and related censoring variable, 
        ## thus G_c(min(T_i,tau)) or G_c(T_i) as weights
        if ((fix.cens.weights==1 & k==0) | (fix.cens.weights==0)) {
            data0 <- data
            time0 <- time
            status0 <- status
        }
        cond0 <- time0>tau
        if (!fix.cens.weights) {
            status0[cond0 & status==1] <- 3 ## Not-censored if T>tau
        }
        data0[,outcome] <- data[outcome]
        data0[cond0,outcome] <- FALSE ## Non-case if T>tau
        if (!fix.cens.weights) time0[cond0] <- tau
        if ((fix.cens.weights & k==0) | (!fix.cens.weights)) {
            
            if (ncol(censtime)==3) { ## truncation...
                suppressWarnings(data0$S <- Surv(ltimes,time0,status0==1))
            } else {
                data0$S <- Surv(time0,status0==1)
            }
            ## data0$status0 <- status0
            ## data0$time0 <- time0
            ## data0$y0 <- data0[,outcome]
            dataw <- ipw(update(cens.formula,S~.), data=subset(data0,ltimes<time0), cens.model=cens.model,
                         weight.name=weights,obs.only=TRUE,cluster=id)            
        }
        if (fix.cens.weights) {
            timevar <- dataw$S[,1]
            dataw[,outcome] <- (dataw[,outcome])*(timevar<tau)
        }
        k <- k+1
        if (ncol(censtime)==3) { ## truncation...
            
        }
        if (return.data) return(dataw)
        if (ncol(censtime)==3) { ## truncation...
            dataw[,weights] <- dataw[,weights]*dataw[,trunc.weights]
        }        
        args <- c(list(x=formula,data=dataw,id=id,weights=weights, pairs.only=pairs.only), list(...))
        suppressWarnings(b <- do.call(estimator, args))
        ## suppressWarnings(b <- biprobit(formula, data=dataw, id=id, weights=weights, pairs.only=pairs.only,...))
        if (length(breaks)>1) res <- c(res,list(summary(b,...)))
    }
    if (length(breaks)==1) {
        return(structure(b,time=breaks))
    }
    if (missing(summary.function)) {
        summary.function <- function(x,...) x$all
    } 
    mycoef <- lapply(rev(res),function(x) summary.function(x))
    res <- list(varname="Time",var=rev(breaks),coef=mycoef,summary=rev(res),call=m,type="time")
    class(res) <- "timemets"
    return(res)    
} ## }}} 


biprobit.time2 <- function(formula,data,id,...,
                          breaks=Inf,pairs.only=TRUE,
                          cens.formula,cens.model="aalen",weights="w") {
## {{{ 

    m <- match.call(expand.dots = FALSE)
    m <- m[match(c("","data"),names(m),nomatch = 0)]
    Terms <- terms(cens.formula,data=data)
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    M <- eval(m,envir=parent.frame())

    censtime <- model.extract(M, "response")
    status <- censtime[,2]
    time <- censtime[,1]
    outcome <- as.character(terms(formula)[[2]])    
    if (is.null(breaks)) breaks <-  quantile(time,c(0.25,0.5,0.75,1))

    outcome0 <- paste(outcome,"_dummy")
    res <- list()
    for (tau in breaks) {
        if (length(breaks)>1) message(tau)
        data0 <- data
        time0 <- time
        cond0 <- time0>tau
        status0 <- status
        status0[cond0 & status==1] <- 3 ## Censored
        data0[cond0,outcome] <- FALSE
        time0[cond0] <- tau
        data0$S <- survival::Surv(time0,status0==1)        
        dataw <- ipw(update(cens.formula,S~.), data=data0, cens.model=cens.model,
                     cluster=id,weight.name=weights,obs.only=TRUE)
        message("control")
        suppressWarnings(b <- biprobit(formula, data=dataw, id=id, weights=weights, pairs.only=pairs.only,...))
        res <- c(res,list(summary(b)))
    }
    if (length(breaks)==1) return(b)
    res <- list(varname="Time",var=breaks,coef=lapply(res,function(x) x$all),summary=res,call=m,type="time")
    class(res) <- "timemets"
    return(res)    
} ## }}} 

biprobit.time.trunc.test <- function(formula,data,id,...,
                          breaks=NULL,n.times=20,pairs.only=TRUE,fix.cens.weights=FALSE,
                          cens.formula,cens.model="aalen",weights="w",messages=FALSE,
                          return.data=FALSE,theta.formula=~1,trunc.weights="w2",
                          estimator="biprobit", summary.function) {
## {{{ 
    
    m <- match.call(expand.dots = FALSE)
    m <- m[match(c("","data"),names(m),nomatch = 0)]
    Terms <- terms(cens.formula,data=data)
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    M <- eval(m,envir=parent.frame())

    censtime <- model.extract(M, "response")
    if (pairs.only) {
        ii <- sort(as.matrix(na.omit(fast.reshape(seq(nrow(data)),id=data[,id]))))
        data <- data[ii,,drop=FALSE]
        censtime <- censtime[ii]
    }

    ltimes <- 0
    if (ncol(censtime)==3) {  ## Calculate probability of not being truncated via Clayton-Oakes (ad hoc combining causes)
        status <- censtime[,3]
        noncens <- !status
        time <- censtime[,2]
        ltimes <- censtime[,1]        
        data$truncsurv <- Surv(ltimes,time,noncens)
        trunc.formula <- update(formula,truncsurv~.)       
        ud.trunc <- aalen(trunc.formula,data=data,robust=0,n.sim=0,residuals=0,silent=1,max.clust=NULL,
                          clusters=data[,id], ...)
        X <- model.matrix(trunc.formula,data)
        ##dependX0 <- model.matrix(theta.formula,data)
        dependX0 <- X       
        twostage.fit <- two.stage(ud.trunc,
                                 data=data,robust=0,detail=0,
                                 theta.des=dependX0)#,Nit=20,step=1.0,notaylor=1)        
        Xnam <- colnames(X)
        ww <- fast.reshape(cbind(X,".num"=seq(nrow(X)),".lefttime"=ltimes),varying=c(".num",".lefttime"),id=data[,id])
        dependX <- as.matrix(ww[,Xnam,drop=FALSE])
        nottruncpair.prob <- predict.two.stage(twostage.fit,X=ww[,Xnam],
                                               times=ww[,".lefttime1"],times2=ww[,".lefttime2"],
                                               theta.des=dependX)$St1t2
        data[,trunc.weights] <- 0
        data[ww[,".num1"],trunc.weights] <- nottruncpair.prob
        data[ww[,".num2"],trunc.weights] <- nottruncpair.prob
        
    } else {
        status <- censtime[,2]
        time <- censtime[,1]
    }

    outcome <- as.character(terms(formula)[[2]])
    jj <- jumptimes(time,data[,outcome],data[,id],sample=n.times)
    lastjump <- tail(jj,1)
    if (is.null(breaks)) {
        breaks <- jj
    }
    if (any(breaks>lastjump)) {
        breaks <- unique(pmin(breaks,tail(jj,1)))        
        if (messages) message("Looking at event before time ",max(breaks))
    }    

    outcome0 <- paste(outcome,"_dummy")
    res <- list(); k <- 0
    breaks <- rev(breaks)
    for (tau in breaks) {
        if (length(breaks)>1 && messages) message(tau)
        ## construct min(T_i,tau) or T_i and related censoring variable, 
        ## thus G_c(min(T_i,tau)) or G_c(T_i) as weights
        if ((fix.cens.weights==1 & k==0) | (fix.cens.weights==0)) {
            data0 <- data
            time0 <- time
            status0 <- status
        }
        cond0 <- time0>tau
        if (!fix.cens.weights) {
            status0[cond0 & status==1] <- 3 ## Not-censored if T>tau
        }
        data0[,outcome] <- data[outcome]
        data0[cond0,outcome] <- FALSE ## Non-case if T>tau
        if (!fix.cens.weights) time0[cond0] <- tau
        if ((fix.cens.weights & k==0) | (!fix.cens.weights)) {
            
            if (ncol(censtime)==3) { ## truncation...
                suppressWarnings(data0$S <- Surv(ltimes,time0,status0==1))
            } else {
                data0$S <- Surv(time0,status0==1)
            }
            ## data0$status0 <- status0
            ## data0$time0 <- time0
            ## data0$y0 <- data0[,outcome]
            dataw <- ipw(update(cens.formula,S~.), data=subset(data0,ltimes<time0), cens.model=cens.model,
                         weight.name=weights,obs.only=TRUE,cluster=id)            
        }
        if (fix.cens.weights) {
            timevar <- dataw$S[,1]
            dataw[,outcome] <- (dataw[,outcome])*(timevar<tau)
        }
        k <- k+1
        if (ncol(censtime)==3) { ## truncation...
            
        }
        if (return.data) return(dataw)
        if (ncol(censtime)==3) { ## truncation...
            dataw[,weights] <- dataw[,weights]*dataw[,trunc.weights]
        }        
        args <- c(list(x=formula,data=dataw,id=id,weights=weights, pairs.only=pairs.only), list(...))
        suppressWarnings(b <- do.call(estimator, args))
        ## suppressWarnings(b <- biprobit(formula, data=dataw, id=id, weights=weights, pairs.only=pairs.only,...))
        if (length(breaks)>1) res <- c(res,list(summary(b,...)))
    }
    if (length(breaks)==1) {
        return(structure(b,time=breaks))
    }
    if (missing(summary.function)) {
        summary.function <- function(x,...) x$all
    } 
    mycoef <- lapply(rev(res),function(x) summary.function(x))
    res <- list(varname="Time",var=rev(breaks),coef=mycoef,summary=rev(res),call=m,type="time")
    class(res) <- "timemets"
    return(res)    
} ## }}} 


