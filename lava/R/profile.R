##' @export
profile.lvmfit <- function(fitted,idx,tau,...) {
  mm <- parfix(Model(fitted),idx,tau)
  index(mm) <- reindex(mm,zeroones=TRUE,deriv=TRUE)
  fixed <- attributes(mm)$fixed
  plogl <- function(tau0) {
    for (i in fixed$v) {
      mm$mean[[i]] <- tau0
    }
    for (i in seq_len(nrow(fixed$A))) {
        index(mm)$A[fixed$A[i,1],fixed$A[i,2]] <-
            mm$fix[fixed$A[i,1],fixed$A[i,2]] <- tau0
    }
    for (i in seq_len(nrow(fixed$P))) {
        index(mm)$P[fixed$P[i,1],fixed$P[i,2]] <-
            mm$covfix[fixed$P[i,1],fixed$P[i,2]] <- tau0
    }
    for (i in length(fixed$e)) {
        index(mm)$exfix[i] <- tau0
    }
    dots <- list(...)
    dots$silent <- FALSE
    if (!is.null(dots$control))
      control <- dots$control
    else
      control <- list()
    control$start <- coef(fitted)
    dots$control <- control
    dots$index <- FALSE
    dots$fix <- FALSE
    dots$silent <- TRUE
    dots$quick <- TRUE
    dots$data <- model.frame(fitted)
    dots$x <- mm
    ee <- do.call("estimate",dots)
    return(logLik(mm,p=ee,data=dots$data))
  }
  val <- sapply(tau,plogl)
  attributes(val) <- NULL
  val
}

profci.lvmfit <- function(x,parm,level=0.95,interval=NULL,curve=FALSE,n=20,lower=TRUE,upper=TRUE,...) {
  ll <- logLik(x)-qchisq(level,1)/2
  pp <- function(tau) (profile.lvmfit(x,parm,tau) - ll)
  tau0 <- coef(x)[parm]
  tau0.sd <- x$vcov[parm,parm]^0.5
  if (is.null(interval)) {
      interval <- tau0 + 3*c(-1,1)*tau0.sd
      if (parm%in%(variances(x)+index(x)$npar.mean))
          interval[1] <- max(1e-5,interval[1])
  }
  if (curve) {
    xx <- seq(interval[1],interval[2],length.out=n)
    val <- sapply(xx,pp)
    res <- cbind(par=xx,val=val)
    return(res)
  }
  low <- up <- NA
  if (lower)
    low <- uniroot(pp,interval=c(interval[1],tau0))$root
  if (upper)
    up <- uniroot(pp,interval=c(tau0,interval[2]))$root
  ##  res <- rbind(lower$root,upper$root); rownames(res) <- coef()
  return(c(low,up))
}
