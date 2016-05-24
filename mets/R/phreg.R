
###{{{ phreg0 

phreg0 <- function(X,entry,exit,status,id=NULL,strata=NULL,beta,stderr=TRUE,method="NR",...) {
  p <- ncol(X)
  if (missing(beta)) beta <- rep(0,p)
  if (p==0) X <- cbind(rep(0,length(exit)))
  if (!is.null(strata)) {
    stratalev <- levels(strata)
    strataidx <- lapply(stratalev,function(x) which(strata==x))
    if (!all(unlist(lapply(strataidx,function(x) length(x)>0))))
      stop("Strata without any observation")
    dd <- lapply(strataidx, function(ii)                      
                 .Call("FastCoxPrep",
                       entry[ii],exit[ii],status[ii],
                       as.matrix(X)[ii,,drop=FALSE],
                       id[ii],
                       !is.null(id),
                       !is.null(entry),
                       package="mets"))
    if (!is.null(id))
      id <- unlist(lapply(dd,function(x) x$id[x$jumps+1]))
    obj <- function(pp,U=FALSE,all=FALSE) {
      val <- lapply(dd,function(d)
                    with(d,
                         .Call("FastCoxPL",pp,X,XX,sign,jumps,package="mets")))
      ploglik <- Reduce("+",lapply(val,function(x) x$ploglik))
      gradient <- Reduce("+",lapply(val,function(x) x$gradient))
      hessian <- Reduce("+",lapply(val,function(x) x$hessian))
      if (all) {
        U <- do.call("rbind",lapply(val,function(x) x$U))
        time <- lapply(dd,function(x) x$time[x$ord+1])
        ord <- lapply(dd,function(x) x$ord+1)
        jumps <- lapply(dd,function(x) x$jumps+1)
        jumptimes <- lapply(dd,function(x) x$time[x$ord+1][x$jumps+1])
        S0 <- lapply(val,function(x) x$S0)
        nevent  <- unlist(lapply(S0,length))
        return(list(ploglik=ploglik,gradient=gradient,hessian=hessian,
                    U=U,S0=S0,nevent=nevent,
                    ord=ord,time=time,jumps=jumps,jumptimes=jumptimes))
      }
      structure(-ploglik,gradient=-gradient,hessian=-hessian)
    }
  } else {
      browser()      
      system.time(dd <- .Call("mets_FastCoxPrep",
                              entry,exit,status,X,
                              as.integer(seq_along(entry)),
                              is.null(id),
                              !is.null(entry),
                              package="mets"))
      if (!is.null(id))
          id <- dd$id[dd$jumps+1]
      obj <- function(pp,U=FALSE,all=FALSE) {
          val <- with(dd,
                      .Call("FastCoxPL",pp,X,XX,sign,jumps,package="mets"))
          if (all) {
              val$time <- dd$time[dd$ord+1]
              val$ord <- dd$ord+1
              val$jumps <- dd$jumps+1
              val$jumptimes <- val$time[val$jumps]
              val$nevent <- length(val$S0)
              return(val)
          }
          with(val, structure(-ploglik,gradient=-gradient,hessian=-hessian))
      }
  }
  opt <- NULL
  if (p>0) {      
      if (tolower(method)=="nr") {
          opt <- lava::NR(beta,obj,...)
          opt$estimate <- opt$par
      } else {
          opt <- nlm(obj,beta,...)
          opt$method <- "nlm"
      }
      cc <- opt$estimate;  names(cc) <- colnames(X)
      if (!stderr) return(cc)
      val <- c(list(coef=cc),obj(opt$estimate,all=TRUE))
  } else {
      val <- obj(0,all=TRUE)
      val[c("ploglik","gradient","hessian","U")] <- NULL
  }
  res <- c(val,
           list(strata=strata,
                entry=entry,
                exit=exit,
                status=status,
                p=p,
                X=X,
                id=id, opt=opt))
  class(res) <- "phreg"
  res
}

###}}} phreg0

###{{{ simcox

##' @export
simCox <- function(n=1000, seed=1, beta=c(1,1), entry=TRUE) {
  if (!is.null(seed))
      set.seed(seed)
  m <- lvm()
  regression(m,T~X1+X2) <- beta
  distribution(m,~T+C) <- coxWeibull.lvm(scale=1/100)
  distribution(m,~entry) <- coxWeibull.lvm(scale=1/10)
  m <- eventTime(m,time~min(T,C=0),"status")
  d <- sim(m,n);
  if (!entry) d$entry <- 0
  else d <- subset(d, time>entry,select=-c(T,C))
  return(d)
}


###}}} simcox

###{{{ phreg
##' Fast Cox PH regression
##'
##' Fast Cox PH regression
##' @param formula formula with 'Surv' outcome (see \code{coxph})
##' @param data data frame
##' @param ... Additional arguments to lower level funtions
##' @author Klaus K. Holst
##' @export
##' @aliases phreg phreg.par
##' @examples
##' simcox <- function(n=1000, seed=1, beta=c(1,1), entry=TRUE) {
##'   if (!is.null(seed))
##'     set.seed(seed)
##'   library(lava)
##'   m <- lvm()
##'   regression(m,T~X1+X2) <- beta
##'   distribution(m,~T+C) <- coxWeibull.lvm(scale=1/100)
##'   distribution(m,~entry) <- coxWeibull.lvm(scale=1/10)
##'   m <- eventTime(m,time~min(T,C=0),"status")
##'   d <- sim(m,n);
##'   if (!entry) d$entry <- 0
##'   else d <- subset(d, time>entry,select=-c(T,C))
##'   return(d)
##' }
##' \dontrun{
##' n <- 1e3;
##' d <- mets:::simCox(n); d$id <- seq(nrow(d)); d$group <- factor(rbinom(nrow(d),1,0.5))
##' 
##' (m1 <- phreg(Surv(entry,time,status)~X1+X2,data=d))
##' (m2 <- coxph(Surv(entry,time,status)~X1+X2+cluster(id),data=d))
##' (coef(m3 <-cox.aalen(Surv(entry,time,status)~prop(X1)+prop(X2),data=d)))
##' 
##' 
##' (m1b <- phreg(Surv(entry,time,status)~X1+X2+strata(group),data=d))
##' (m2b <- coxph(Surv(entry,time,status)~X1+X2+cluster(id)+strata(group),data=d))
##' (coef(m3b <-cox.aalen(Surv(entry,time,status)~-1+group+prop(X1)+prop(X2),data=d)))
##' 
##' m <- phreg(Surv(entry,time,status)~X1*X2+strata(group)+cluster(id),data=d)
##' m
##' plot(m,ylim=c(0,1))
##' }
phreg <- function(formula,data,...) {
  cl <- match.call()
  m <- match.call(expand.dots = TRUE)[1:3]
  special <- c("strata", "cluster")
  Terms <- terms(formula, special, data = data)
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  Y <- model.extract(m, "response")
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
  id <- strata <- NULL
  if (!is.null(attributes(Terms)$specials$cluster)) {
    ts <- survival::untangle.specials(Terms, "cluster")
    Terms  <- Terms[-ts$terms]
    id <- m[[ts$vars]]
  }
  if (!is.null(stratapos <- attributes(Terms)$specials$strata)) {
    ts <- survival::untangle.specials(Terms, "strata")
    Terms  <- Terms[-ts$terms]
    strata <- m[[ts$vars]]
  }  
  X <- model.matrix(Terms, m)
  if (!is.null(intpos  <- attributes(Terms)$intercept))
    X <- X[,-intpos,drop=FALSE]
  if (ncol(X)==0) X <- matrix(nrow=0,ncol=0)

  res <- c(phreg0(X,entry,exit,status,id,strata,...),list(call=cl,model.frame=m))
  class(res) <- "phreg"
  
  res
}
###}}} phreg

###{{{ vcov

##' @export
vcov.phreg  <- function(object,...) {    
  res <- crossprod(ii <- iid(object,...))
  attributes(res)$ncluster <- attributes(ii)$ncluster
  attributes(res)$invhess <- attributes(ii)$invhess
  colnames(res) <- rownames(res) <- names(coef(object))
  res
}

###}}} vcov

###{{{ coef

##' @export
coef.phreg  <- function(object,...) {
  object$coef
}

###}}} coef

###{{{ iid

##' @export
iid.phreg  <- function(x,...) {
    invhess <- solve(x$hessian)
  ncluster <- NULL
  if (!is.null(x$id)) {
    ii <- mets::cluster.index(x$id)
    UU <- matrix(nrow=ii$uniqueclust,ncol=ncol(invhess))
    for (i in seq(ii$uniqueclust)) {
      UU[i,] <- colSums(x$U[ii$idclustmat[i,seq(ii$cluster.size[i])]+1,,drop=FALSE])
    }
    ncluster <- nrow(UU)
  } else {
      UU <- x$U
  }
  structure(UU%*%invhess,invhess=invhess,ncluster=ncluster)
}

###}}}

###{{{ summary

##' @export
summary.phreg <- function(object,se="robust",...) {
  cc <- ncluster <- NULL
  if (length(object$p)>0) {
    I <- -solve(object$hessian)
    V <- vcov(object)
    cc <- cbind(coef(object),diag(V)^0.5,diag(I)^0.5)
    cc  <- cbind(cc,2*(pnorm(abs(cc[,1]/cc[,2]),lower.tail=FALSE)))
    colnames(cc) <- c("Estimate","S.E.","dU^-1/2","P-value")
    if (!is.null(ncluster <- attributes(V)$ncluster))
    rownames(cc) <- names(coef(object))
  }
  Strata <- levels(object$strata)
  if (!is.null(Strata)) {
    n <- unlist(lapply(object$time,length))
  } else {
    n <- length(object$time)    
  }  
  res <- list(coef=cc,n=n,nevent=object$nevent,
              strata=Strata,ncluster=ncluster)
  class(res) <- "summary.phreg"
  res
}

###}}} summary

###{{{ print.summary

##' @export
print.summary.phreg  <- function(x,max.strata=5,...) {
  cat("\n")
  nn <- cbind(x$n, x$nevent)
  rownames(nn) <- levels(x$strata); colnames(nn) <- c("n","events")
  if (is.null(rownames(nn))) rownames(nn) <- rep("",NROW(nn))
  if (length(x$strata)>max.strata) {
      nn <- rbind(c(colSums(nn),length(x$strata)));
      colnames(nn) <- c("n","events","stratas")
      rownames(nn) <- ""
  } 
  print(nn,quote=FALSE)  
  if (!is.null(x$ncluster)) cat("\n ", x$ncluster, " clusters\n",sep="")
  if (!is.null(x$coef)) {
    cat("\n")
    printCoefmat(x$coef,...)
  }
  cat("\n")
}

###}}} print.summary

###{{{ predict

predictPhreg <- function(jumptimes,S0,beta,time=NULL,X=NULL,surv=FALSE,...) {
    ## Brewslow estimator
    chaz <- cbind(jumptimes,cumsum(1/S0))
    if (!is.null(time)) {
        chaz <- Cpred(chaz,time)
    }
    colnames(chaz) <- c("time","chaz")
    if (!is.null(X)) {
      H <- exp(X%*%beta)
      if (nrow(chaz)==length(H)) {
        chaz[,2] <- chaz[,2]*H
      } else {
        chaz2 <- c()
        X <- rbind(X)
        for (i in seq(nrow(X)))
          chaz2 <- rbind(chaz2,
                         cbind(chaz[,1],chaz[,2]*H[i],
                               rep(1,nrow(chaz))%x%X[i,,drop=FALSE]))
        chaz <- chaz2;
        nn <- c("time","chaz",names(beta))
        colnames(chaz) <- nn
      }
    }
    if (surv) {    
      chaz[,2] <- exp(-chaz[,2])
      colnames(chaz)[2] <- "surv"
    }
    return(chaz)
}

##' @export
predict.phreg  <- function(object,data,surv=FALSE,time=object$exit,X=object$X,strata=object$strata,...) {
    if (!is.null(object$strata)) {
        lev <- levels(object$strata)
        if (!is.null(object$strata) &&
            !(is.list(time) & !is.data.frame(time)) &&
            !(is.list(X) & !is.data.frame(X))) {
            X0 <- X
            time0 <- time
            X <- time <- c()
            for (i in seq(length(lev))) {
                idx <- which(strata==lev[i])
                X <- c(X,list(X0[idx,,drop=FALSE]))
                time <- c(time,list(time0[idx]))
            }
        }
        chaz <- c()
        for (i in seq(length(lev)))
            chaz <- c(chaz,list(predictPhreg(object$jumptimes[[i]],
                                             object$S0[[i]],
                                             coef(object),
                                             time[[i]],X[[i]],surv)))
        names(chaz) <- lev    
    } else {
        chaz <- predictPhreg(object$jumptimes,object$S0,coef(object),time,X,surv)
    }
    return(chaz)
}

###}}} predict

###{{{ plot

##' @export
plot.phreg  <- function(x,surv=TRUE,X=NULL,time=NULL,add=FALSE,...) {
    if (!is.null(X) && nrow(X)>1) {
        P <- lapply(split(X,seq(nrow(X))),function(xx) predict(x,X=xx,time=time,surv=surv))
    } else {
        P <- predict(x,X=X,time=time,surv=surv)
    }
    if (!is.list(P)) {
        if (add) {
            lines(P,type="s",...)
        } else {
            plot(P,type="s",...)
        }
        return(invisible(P))
    }

    if (add) {
        lines(P[[1]][,1:2],type="s",lty=1,col=1,...)
    } else {
        plot((P[[1]])[,1:2],type="s",lty=1,col=1,...)
    }
    for (i in seq_len(length(P)-1)+1) {
        lines(P[[i]][,1:2],type="s",lty=i,col=i,...)   
    }
    return(invisible(P))
}

###}}} plot

###{{{ print
##' @export
print.phreg  <- function(x,...) {
  cat("Call:\n")
  dput(x$call)
  print(summary(x),...)
}
###}}} print


