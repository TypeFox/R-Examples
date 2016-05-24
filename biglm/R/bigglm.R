bigglm<-function(formula, data, family=gaussian(),...)
    UseMethod("bigglm", data)
setGeneric("bigglm", signature=c("formula","data"))


bigglm.data.frame<-function(formula, data, ..., chunksize=5000){
    n<-nrow(data)
    cursor<-0
    datafun<-function(reset=FALSE){
        if (reset){
            cursor<<-0
            return(NULL)
        }
        if (cursor>=n)
            return(NULL)
        start<-cursor+1
        cursor<<-cursor+min(chunksize, n-cursor)
        data[start:cursor,]
    }
    rval<-bigglm(formula=formula, data=datafun,...)
    rval$call<-sys.call()
    rval$call[[1]]<-as.name(.Generic)
    rval
}

bigglm.function<-function(formula, data, family=gaussian(), weights=NULL,
                            sandwich=FALSE, maxit=8, tolerance=1e-7,
                            start=NULL, quiet=FALSE,...){

    tt<-terms(formula)
    beta <- start
    etafun <- function(x) if(is.null(beta)) rep(0,nrow(x)) else x%*%beta
    
    
    converged<-FALSE
    for (i in 1:maxit){
        firstchunk <- TRUE
        deviance<-0
        rss<-0
        data(reset=TRUE)
        n<-0
        while(!is.null(chunk<-data(reset=FALSE))){
            n<-n+nrow(chunk)
            mf<-model.frame(tt,chunk)
            mm<-model.matrix(tt,mf)
            p<-NCOL(mm)
            if (!is.null(weights)){
                if (!inherits(weights, "formula"))
                    stop("`weights' must be a formula")
                w<-model.frame(weights, chunk)[[1]]
            } else w<-rep(1,nrow(mm))
            if (firstchunk) {
                qr<-bigqr.init(p)
                assn<-attr(mm,"assign")
                if(sandwich)
                    xyqr<-bigqr.init(p*(p+1))
            }
            if (!identical(assn, attr(mm,"assign")))
                stop("model matrices incompatible")
            y<-model.response(mf)
            if(is.null(off<-model.offset(mf))) off<-0
            eta<-etafun(mm)+off
            mu <- family$linkinv(eta)
            dmu <- family$mu.eta(eta)
            z<- eta+(y-mu)/dmu 
            ww<-w*dmu*dmu/(family$variance(mu))
            qr<-update(qr,mm,z-off,ww)
            if(!is.null(beta)){
                deviance<-deviance+sum(family$dev.resids(y,mu,w))
                rss<-rss+sum((y-mu)^2/(w*family$variance(mu)))*(sum(w)/length(w))
                if (sandwich){
                    xx<-matrix(nrow=nrow(mm), ncol=p*(p+1))
                    xx[,1:p]<-mm*(drop(z)-off)
                    for(i in 1:p)
                        xx[,p*i+(1:p)]<-mm*mm[,i]
                    xyqr<-update(xyqr,xx,rep(0,nrow(mm)),ww*ww)
                }
            }
            firstchunk <- FALSE
        }
        iwlm <- list(call=sys.call(-1), qr=qr, iterations=i,
                   assign=attr(mm,"assign"), terms=tt, converged=FALSE,
                   n=n,names=colnames(mm), weights=weights,rss=rss)
        if(sandwich)
            iwlm$sandwich <- list(xy=xyqr)
        class(iwlm) <- "biglm"

        betaold <- beta
        beta <- coef(iwlm)
        if (i >= maxit){
          if (!quiet) warning("ran out of iterations and failed to converge")
            break
          }
            
        if (!is.null(betaold)){
            delta <- (betaold-beta)/sqrt(diag(vcov(iwlm)))
            if (max(abs(delta)) < tolerance){
                iwlm$converged<-TRUE
                break
            }
        }

        
    }
    
    rval <- iwlm
    rval$family <- family
    rval$deviance <- deviance
    rval$df.resid <- rval$n-length(rval$qr$D)
    class(rval) <- c("bigglm","biglm")
    rval
}




print.bigglm<-function(x,...){
  cat("Large data regression model: ")
  print(x$call)
  cat("Sample size = ",x$n,"\n")
  if (is.null(x$converged)|| !x$converged)
    cat("failed to converge after", x$iterations," iterations\n")
  invisible(x)
}

