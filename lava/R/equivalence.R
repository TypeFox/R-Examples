##' Identify candidates of equivalent models
##'
##' Identifies candidates of equivalent models
##'
##'
##' @param x \code{lvmfit}-object
##' @param rel Formula or character-vector specifying two variables to omit from
##' the model and subsequently search for possible equivalent models
##' @param tol Define two models as empirical equivalent if the absolute
##' difference in score test is less than \code{tol}
##' @param k Number of parameters to test simultaneously. For \code{equivalence}
##' the number of additional associations to be added instead of \code{rel}.
##' @param omitrel if \code{k} greater than 1, this boolean defines wether to
##' omit candidates containing \code{rel} from the output
##' @param \dots Additional arguments to be passed to the low level functions
##' @author Klaus K. Holst
##' @seealso \code{\link{compare}}, \code{\link{modelsearch}}
##' @export
equivalence <- function(x,rel,tol=1e-3,k=1,omitrel=TRUE,...) {
  if (missing(rel)) stop("Specify association 'rel' (formula or character vector)")
  if (inherits(rel,"formula")) {
    myvars <- all.vars(rel)
  } else {
    myvars <- rel
  }
  if (length(myvars)!=2) stop("Two variables only")
  x0 <- Model(x)
  cancel(x0) <- rel
  e0 <- estimate(x0,data=model.frame(x),weight=Weight(x),estimator=x$estimator,...)
  if (k!=1) {
    p0 <- coef(x)
    p0[] <- 0
    p0[match(names(coef(e0)),names(p0))] <- coef(e0)
    S0 <- score(x,p=p0)[,,drop=TRUE];
    I0 <- information(x,p=p0)
    T0 <- rbind(S0)%*%solve(I0)%*%cbind(S0); names(T0) <- "Q"
  }
  s <- modelsearch(e0,k=k,...)
  relname <- c(paste(myvars,collapse=lava.options()$symbol[2]),
               paste(rev(myvars),collapse=lava.options()$symbol[2]))
  relidx <- NULL
  if (k==1) {
    relidx <- na.omit(match(relname,s$res[,"Index"]))
    T0 <- s$test[relidx,1]
  }
  T <- s$test[,1]
  Equiv <- setdiff(which(abs(T-T0)<tol),relidx)
  Improve <- which((T-T0)>tol)
  if (omitrel) { ## Don't save models including 'rel'
    keep <- c()
    if (length(Equiv)>0) {
      for (i in seq_len(length(Equiv))) {
        newvars <- s$var[[Equiv[i]]]
        if (!any(apply(newvars,1,function(z) all(z%in%myvars)))) keep <- c(keep,Equiv[i])
      }
      Equiv <- keep
    }
    keep <- c()
    if (length(Improve)>0) {
      for (i in seq_len(length(Improve))) {
        newvars <- s$var[[Improve[i]]]
        if (!any(apply(newvars,1,function(z) all(z%in%myvars)))) keep <- c(keep,Improve[i])
      }
      Improve <- keep
    }
  }
  eqvar <- ivar <- NULL
  models <- list()
  if (length(Equiv)>0){
    for (i in seq_len(length(Equiv))) {
      xnew <- x0
      newvars <- s$var[[Equiv[i]]]
      for (j in seq_len(nrow(newvars))) {
        exo.idx <- which(newvars[j,]%in%index(x0)$exogenous)
        if (length(exo.idx)>0) {
          xnew <- regression(xnew,from=newvars[j,exo.idx],to=newvars[j,setdiff(1:2,exo.idx)])
        } else {
          covariance(xnew) <- newvars
        }
      }
      models <- c(models,list(xnew))
    }
    eqvar <- s$var[Equiv]
  }
  if (length(Improve)>0)   {
      for (i in seq_len(length(Improve))) {
      xnew <- x0
      newvars <- s$var[[Improve[i]]]
      for (j in seq_len(nrow(newvars))) {
        exo.idx <- which(newvars[j,]%in%index(x0)$exogenous)
        if (length(exo.idx)>0) {
          xnew <- regression(xnew,from=newvars[j,exo.idx],to=newvars[j,setdiff(1:2,exo.idx)])
        } else {
          covariance(xnew) <- newvars
        }
      }
      models <- c(models,list(xnew))
    }
    ivar <- s$var[Improve]
  }
  res <- list(equiv=eqvar, improve=ivar, scoretest=s, models=models, I=Improve, E=Equiv, T0=T0, vars=myvars)
  class(res) <- "equivalence"
  return(res)
}

##' @export
print.equivalence <- function(x,...) {
  cat("  0)\t ",paste0(x$vars,collapse=lava.options()$symbol[2]),"  (",formatC(x$T0),")\n")
  cat("Empirical equivalent models:\n")
  if (length(x$E)==0)
    cat("\t none\n")
  else
    for (i in seq_len(length(x$E))) {
      cat("  ",i,")\t ",  x$scoretest$res[x$E[i],"Index"],
          "  (",x$scoretest$res[x$E[i],1],")",
          "\n",sep="")
    }
  cat("Candidates for model improvement:\n")
  if (length(x$I)==0)
    cat("\t none\n")
  else
  for (i in seq_len(length(x$I))) {
      cat("  ",i,")\t ",  x$scoretest$res[x$I[i],"Index"],
          "  (",x$scoretest$res[x$I[i],1],")",
          "\n",sep="")
  }
  invisible(x)
}

holm <- function(p) {
  k <- length(p)
  w <- 1/k
  ii <- order(p)
  po <- p[ii]
  qs <- min(1,po[1]/w)
  for (i in 2:k) {
      qs <- c(qs, min(1, max(qs[i-1],po[i]*(1-w*(i-1))/w)))
    }
  return(qs)
}
