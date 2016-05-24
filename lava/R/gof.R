##' @export
rsq <- function(x,stderr=FALSE) {

    if (stderr) {
        v <- endogenous(x)
        vpar <- paste(v,v,sep=lava.options()$symbol[2])
        iid.v <- iid(model.frame(x)[,v])
        iid.mod <- iid(x)
        coef0 <- c(attributes(iid.v)$coef[vpar],
                   coef(x)[vpar])
        iid0 <- cbind(iid.v[,vpar],iid.mod[,vpar])
        p <- length(v)
        idx <- seq_len(p);
        ee <- estimate(NULL,data=NULL,
                       function(x) {
                           res <- (x[idx]-x[idx+p])/x[idx]
                           names(res) <- v
                           as.list(res)
                       },
                       print=function(x,...) {
                           cat("\nR-squared:\n\n")
                           print(x$coefmat)
                       },
                       coef=coef0, iid=iid0)

        res <- ee
        ##res <- list(ee)
        ## for (lat in latent(x)) {

        ##     v <- intersect(children(x,lat),endogenous(x))
        ##     vpar <- paste(v,v,sep=lava.options()$symbol[2])
        ##     lpar <- paste(lat,lat,sep=lava.options()$symbol[2])
        ##     rpar <- paste(v,lat,sep=lava.options()$symbol[1])
        ##     fix <- c(x$model$fix[lat,v,drop=TRUE],x$model$covfix[lat,lat])
        ##     pp <- coef(x)
        ##     idx <- x$model$parpos$A[lat,v]
        ##     idx2 <- x$model$parpos$P[lat,lat]
        ##     p0 <- c(idx,idx2)
        ##     p1 <- setdiff(unique(p0),0)
        ##     p2 <- match(p0,p1)

        ##     k <- length(v)
        ##     coef0 <- c(pp[p1],attributes(iid.v)$coef[vpar])
        ##     iid0 <- cbind(iid.mod[,p1],iid.v[,vpar])
        ##     ee <- estimate(NULL,data=NULL,
        ##                    function(p) {
        ##                        p. <- p[p2]
        ##                        p.[is.na(p.)] <- fix[is.na(p.)]
        ##                        res <- p.[seq_len(k)]^2*p.[k+1]/tail(p,k)
        ##                        names(res) <- v
        ##                        as.list(res)
        ##                    },
        ##                    print=function(x,...) {
        ##                        cat("\nVariance explained by '", lat,"':\n\n",sep="")
        ##                        print(x$coefmat)
        ##                    },coef=coef0,iid=iid0)
        ##     res <- c(res,list(ee))
        ## }

        return(res)
    }



    v <- c(endogenous(x),setdiff(latent(x),parameter(Model(x))))
    res <- coef(x,9,std="yx")
    idx <- with(attributes(res),
                which(type=="variance" & (var==from)))
    nam <- attributes(res)$var[idx]
    res <- 1-res[idx,5]
    names(res) <- nam
    res <- list("R-squared"=res)
    ## M <- moments(x,coef(x))
    ## v <- setdiff(vars(x),exogenous(x))
    ## vvar <- M$Cfull[cbind(v,v)]
    ## rsq <- (vvar-M$P[cbind(v,v)])/vvar

    if (length(latent(x))>0) {
        M <- moments(x,coef(x))
        nn <- names(res)
        for (lat in latent(x)) {
            v <- intersect(children(x,lat),endogenous(x))
            varl <- M$Cfull[lat,lat]
            varv <- M$Cfull[cbind(v,v)]
            rpar <- paste(v,lat,sep=lava.options()$symbol[1])
            fix <- c(x$model$fix[lat,v,drop=TRUE])
            pp <- coef(x)
            idx1 <- x$model$parpos$A[lat,v]
            ##idx2 <- x$model$parpos$P[lat,lat]
            ##idx3 <- x$model$parpos$P[cbind(v,v)]
            p0 <- c(idx1)
            p1 <- setdiff(unique(p0),0)
            p2 <- match(p0,p1)
            p <- coef(x)[p1]
            p. <- p[p2]
            p.[is.na(p.)] <- fix[is.na(p.)]
            k <- length(v)
            val <- (p.^2*varl)/varv; names(val) <- v
            res <- c(res,list(val))
            nn <- c(nn,paste0("Variance explained by '",lat,"'"))
        }
        names(res) <- nn
    }
    res
}

satmodel <- function(object,logLik=TRUE,data=model.frame(object),
                     control=list(trace=1),
                     weight=Weight(object),estimator=object$estimator,
                     missing=inherits(object,"lvm.missing"),
                     regr=FALSE,
                     ...) {
  if (object$estimator=="gaussian" & logLik & !missing) {
    if (class(object)[1]%in%c("multigroupfit","multigroup")) {

      ll <- structure(0,nall=0,nobs=0,df=0,class="logLik")
      for (i in seq_len(Model(object)$ngroup)) {
        l0 <- logLik(Model(Model(object))[[i]],data=model.frame(object)[[i]],type="sat")

        ll <- ll+l0
        for (atr in c("nall","nobs","df"))
          attributes(ll)[[atr]] <- attributes(ll)[[atr]]+attributes(l0)[[atr]]
      }

    }
    return(logLik(object, type="sat"))
  }
  covar <- exogenous(object)
  y <- endogenous(object)
  m0 <- Model(object)
  if (length(covar)>0)
    suppressWarnings(m0 <- regression(m0,y,covar))
  if (length(latent(m0))>0)
    kill(m0) <- latent(m0)
  cancel(m0) <- y
  if (!regr)
    suppressWarnings(covariance(m0) <- y)
  else {
    if (length(y)>1) {
      for (i in seq_len(length(y)-1))
       for (j in seq(i+1,length(y))) {
         m0 <- regression(m0,y[i],y[j])
       }
    }
    exogenous(m0) <- covar
  }
  if (is.null(control$start)) {
    mystart <- rep(0,with(index(m0), npar.mean+npar))
    mystart[variances(m0,mean=TRUE)] <- 1
    control$start <- mystart
  }
  message("Calculating MLE of saturated model:\n")
  e0 <- estimate(m0,data=data,weight=weight,estimator=estimator,silent=TRUE,control=control,missing=missing,...)
  if (logLik)
    return(logLik(e0))
  return(e0)
}

condition <- function(A) {
  suppressWarnings(with(eigen(A),tail(values,1)/head(values,1)))
}



##' Extract model summaries and GOF statistics for model object
##'
##' Calculates various GOF statistics for model object including global
##' chi-squared test statistic and AIC. Extract model-specific mean and variance
##' structure, residuals and various predicitions.
##'
##'
##' @aliases gof gof.lvmfit moments moments.lvm information information.lvmfit
##' score score.lvmfit logLik.lvmfit
##' @param object Model object
##' @param x Model object
##' @param p Parameter vector used to calculate statistics
##' @param data Data.frame to use
##' @param weight2 Optional second data.frame (only for censored observations)
##' @param weight Optional weight matrix
##' @param n Number of observations
##' @param conditional If TRUE the conditional moments given the covariates are
##' calculated. Otherwise the joint moments are calculated
##' @param model String defining estimator, e.g. "gaussian" (see
##' \code{estimate})
##' @param debug Debugging only
##' @param chisq Boolean indicating whether to calculate chi-squared
##' goodness-of-fit (always TRUE for estimator='gaussian')
##' @param level Level of confidence limits for RMSEA
##' @param rmsea.threshold Which probability to calculate, Pr(RMSEA<rmsea.treshold)
##' @param all Calculate all (ad hoc) FIT indices: TLI, CFI, NFI, SRMR, ...
##' @param \dots Additional arguments to be passed to the low level functions
##' @usage
##'
##' gof(object, ...)
##'
##' \method{gof}{lvmfit}(object, chisq=FALSE, level=0.90, rmsea.threshold=0.05,all=FALSE,...)
##'
##' moments(x,...)
##'
##' \method{moments}{lvm}(x, p, debug=FALSE, conditional=FALSE, data=NULL, ...)
##'
##' \method{logLik}{lvmfit}(object, p=coef(object),
##'                       data=model.frame(object),
##'                       model=object$estimator,
##'                       weight=Weight(object),
##'                       weight2=object$data$weight2,
##'                           ...)
##'
##' \method{score}{lvmfit}(x, data=model.frame(x), p=pars(x), model=x$estimator,
##'                    weight=Weight(x), weight2=x$data$weight2, ...)
##'
##' \method{information}{lvmfit}(x,p=pars(x),n=x$data$n,data=model.frame(x),
##'                    model=x$estimator,weight=Weight(x), weight2=x$data$weight2, ...)
##'
##' @return A \code{htest}-object.
##' @author Klaus K. Holst
##' @keywords methods models
##' @export
##' @examples
##' m <- lvm(list(y~v1+v2+v3+v4,c(v1,v2,v3,v4)~x))
##' set.seed(1)
##' dd <- sim(m,1000)
##' e <- estimate(m, dd)
##' gof(e,all=TRUE,rmsea.threshold=0.05,level=0.9)
##'
##'
##' set.seed(1)
##' m <- lvm(list(c(y1,y2,y3)~u,y1~x)); latent(m) <- ~u
##' regression(m,c(y2,y3)~u) <- "b"
##' d <- sim(m,1000)
##' e <- estimate(m,d)
##' rsq(e)
##' ##'
##' rr <- rsq(e,TRUE)
##' rr
##' estimate(rr,contrast=rbind(c(1,-1,0),c(1,0,-1),c(0,1,-1)))
##'
`gof` <-
  function(object,...) UseMethod("gof")

##' @export
gof.lvmfit <- function(object,chisq=FALSE,level=0.90,rmsea.threshold=0.05,all=FALSE,...) {
  n <- object$data$n
  if (class(object)[1]=="multigroupfit") n <- sum(unlist(lapply(object$model$data,nrow)))
  loglik <- logLik(object,...)

  df <- attributes(loglik)$df
  nobs <- attributes(loglik)$nall*length(endogenous(object))
  myAIC <- -2*(loglik - df); attributes(myAIC) <- NULL
  myBIC <- -2*loglik + df*log(nobs); attributes(myBIC) <- NULL

  xconstrain <- intersect(unlist(lapply(constrain(object),function(z) attributes(z)$args)),manifest(object))

  l2D <- mean(object$opt$grad^2)
  S <- vcov(object)
  rnkV <- tryCatch(qr(S)$rank,error=function(...) 0)
  minSV <- attr(S,"minSV")
  condnum <- tryCatch(condition(vcov(object)),error=function(...) NULL)

##  if (class(object)[1]=="lvmfit" & (object$estimator=="gaussian" | chisq) & length(xconstrain)==0 ) {
  if (((object$estimator=="gaussian" & class(object)[1]!="lvm.missing") | chisq) & length(xconstrain)==0 ) {
    res <- list(fit=compare(object), n=n, logLik=loglik, BIC=myBIC, AIC=myAIC)
    q <- res$fit$statistic
    qdf <- res$fit$parameter
    if (all) {
      m0 <- lvm(manifest(object)); exogenous(m0) <- NULL
      e0 <- estimate(m0,model.frame(object))
      g0 <- gof(e0)
      logLikbaseline <- g0$logLik
      qbaseline <- g0$fit$statistic
      qdfbaseline <- g0$fit$parameter
      CFI <- ((qbaseline-qdfbaseline) - (q-qdf))/(qbaseline-qdfbaseline)
      NFI <- (qbaseline-q)/qbaseline
      TLI <- (qbaseline/qdfbaseline-q/qdf)/(qbaseline/qdfbaseline-1)
      S <- object$data$S
      mu <- object$data$mu
      C <- modelVar(object)$C
      xi <- as.vector(modelVar(object)$xi)
      if (is.null(S)) S <- cov(model.frame(object))
      if (is.null(mu)) mu <- colMeans(model.frame(object))
      L <- diag(S)^0.5
      idx <- index(object)$endo.idx
      R <- (diag(1/L))%*%(S-C)%*%(diag(1/L))
      R2 <- (mu-xi)/L
      SRMR <- mean(c(R[upper.tri(R,diag=TRUE)],R2)^2)^0.5
      res <- c(res,list(CFI=CFI,NFI=NFI,TLI=TLI,C=C,S=S,SRMR=SRMR))
      ## if (length(latent(object))>0) {
      ##   SRMR.endo <- mean(c(R[idx,idx][upper.tri(R[idx,idx],diag=TRUE)],R2[idx])^2)^0.5
      ##   res <- c(res,list("SRMR(endogenous)"=SRMR.endo))
      ## }
    }
    ##    if (class(object)[1]=="lvmfit")
    if (rnkV==ncol(vcov(object)) && (!is.null(minSV) && minSV>1e-12)) {

      rmseafun <- function(...) {
        epsilon <- function(lambda) sapply(lambda,function(x)
                                           ifelse(x>0 & qdf>0,sqrt(x/(qdf*(n))),0)) ## n-1,n vs. n-df
        opf <- function(l,p) suppressWarnings(p-pchisq(q,df=qdf,ncp=l))
        ## pchisq(... lower.tail=FALSE)-1
        alpha <- (1-level)/2
        RMSEA <- epsilon(q-qdf)
        B <- max(q-qdf,0)
        lo <- hi <- list(root=0)
        if (RMSEA>0 && opf(0,p=1-alpha)<0) {
          hi <- uniroot(function(x) opf(x,p=1-alpha),c(0,B))
        }
        if (opf(B,p=alpha)<0) {
          lo <- uniroot(function(x) opf(x,p=alpha),c(B,n))
        }
        ci <- c(epsilon(c(hi$root,lo$root)))
        RMSEA <- c(RMSEA=RMSEA,ci);
        names(RMSEA) <- c("RMSEA",paste0(100*c(alpha,(1-alpha)),"%"))
        pval <- pchisq(q,qdf,(n*qdf*rmsea.threshold^2),lower.tail=FALSE)
        res <- list(aa=((q-qdf)/(2*qdf)^0.5),RMSEA=RMSEA, level=level, rmsea.threshold=rmsea.threshold, pval.rmsea=pval)
        return(res)
      }
      rmseaval <- tryCatch(rmseafun(),error=function(e) NULL)
      res <- c(res,rmseaval)
    }
  } else {
    res <- list(n=n, logLik=loglik, BIC=myBIC, AIC=myAIC)
  }

  res <- c(res, L2score=l2D, rankV=rnkV, cond=condnum, k=nrow(vcov(object)))
  class(res) <- "gof.lvmfit"
  return(res)
}

##' @export
print.gof.lvmfit <- function(x,optim=TRUE,...) {
  if (!is.null(x$n)) {
    with(x,
         cat("\n Number of observations =", n, "\n"))
  }
  if (is.null(x$fit)) {
    with(x,
         cat(" Log-Likelihood =", logLik, "\n"))
  }
  with(x,  cat(" BIC =", BIC, "\n",
               "AIC =", AIC, "\n"))
  if (!is.null(x$fit))
  with(x,
       cat(" log-Likelihood of model =", fit$estimate[1], "\n\n",
           "log-Likelihood of saturated model =", fit$estimate[2], "\n",
           "Chi-squared statistic: q =", fit$statistic,
           ", df =", fit$parameter,
           "\n  P(Q>q) =", fit$p.value, "\n"))
  if (!is.null(x$RMSEA)) {
    rr <- round(x$RMSEA*10000)/10000
      rmsea <- paste0(rr[1]," (",rr[2],";",rr[3],")")
    cat("\n RMSEA (",x$level*100,"% CI): ", rmsea,"\n",sep="")
    cat("  P(RMSEA<",x$rmsea.threshold,")=",  x$pval.rmsea,"\n",sep="")
  }
  for (i in c("TLI","CFI","NFI","SRMR","SRMR(endogenous)"))
    if (!is.null(x[[i]])) cat("", i,"=",x[[i]],"\n")

  if (optim) {
    cat("\nrank(Information) = ",x$rankV," (p=", x$k,")\n",sep="")
    cat("condition(Information) = ",x$cond,"\n",sep="")
    cat("mean(score^2) =",x$L2score,"\n")
  }

  invisible(x)
}



## gof.multigroupfit <- function(object,...) {
##   L0 <- logLik(object); df0 <- attributes(L0)$df
##   L1 <- logLik(object,type="sat"); df1 <- attributes(L1)$df

##   df <- df1-df0; names(df) <- "df"
##   Q <- -2*(L0-L1); attributes(Q) <- NULL; names(Q) <- "chisq";
##   pQ <- pchisq(Q,df,lower.tail=FALSE)
##   values <- c(L0,L1); names(values) <- c("log likelihood (model)", "log likelihood (saturated model)")
##   res <- list(statistic = Q, parameter = df,
##               p.value=pQ, method = "Likelihood ratio test",
##               estimate = values)
##   class(res) <- "htest"
##   return(res)
## }
