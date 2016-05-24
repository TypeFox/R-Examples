orlmSet <- function(formula, data, set, direction="increase", n=NULL, base=1, control=orlmcontrol()){
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  x <- model.matrix(mt, mf, contrasts)
  d <- ncol(x)
  if (is.null(n)){
    n <- rep(1,d)
    names(n) <- colnames(x)
  }
  if (length(n) != d) stop("n has not the same length as there are columns in the design matrix!")
  if (is.character(set)) set <- constrSet(n, set=set, direction=direction, base=base)
  out <- lapply(set, function(s){
    orlm(formula, data, constr=s$constr, rhs=s$rhs, nec=s$nec, control=control)
  })
  class(out) <- c("list","orlmlist")
  return(out)
}


orglsSet <- function(formula, data, weights=NULL, correlation=NULL, set, direction="increase", n=NULL, base=1, control=orlmcontrol()){
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  x <- model.matrix(mt, mf, contrasts)
  d <- ncol(x)
  if (is.null(n)){
    n <- rep(1,d)
    names(n) <- colnames(x)
  }
  if (length(n) != d) stop("n has not the same length as there are columns in the design matrix!")
  if (is.character(set)) set <- constrSet(n, set=set, direction=direction, base=base)
  out <- lapply(set, function(s){
    orgls(formula, data, constr=s$constr, rhs=s$rhs, nec=s$nec, control=control, weights=weights, correlation=correlation)
  })
  class(out) <- c("list","orglslist")
  return(out)
}


constrMat <- function(n, type=c("monotone","control","average","laverage","uaverage","caverage"), base=1){
  d <- length(n)
  if (d < 2) stop("less than two groups")
  if (!is.numeric(n)) stop(sQuote("n"), " is not numeric")
  if (base < 1 || base > d) stop("base is not between 1 and ", d)
  if (!is.null(names(n))) varnames <- names(n) else varnames <- 1:d
  type <- match.arg(type)
  switch(type, monotone = {
    cm <- cbind(-diag(d-1),0) + cbind(0,diag(d-1))
  }, control = {
    cm <- diag(d)[-base,]
    cm[,base] <- -1
  }, average = {
    cm0 <- matrix(nrow=d-1, ncol=d)
    cm <- t(sapply(1:(d-1), function(i){
      cm0[i, 1:i] <- -n[1:i]/sum(n[1:i])
      cm0[i, (i+1):d] <- n[(i+1):d]/sum(n[(i+1):d])
      return(cm0[i,])
    }))
  }, laverage = {
    cm0 <- cbind(0,diag(d-1))
    cm <- t(sapply(1:(d-1), function(i){
      cm0[i, 1:i] <- -n[1:i]/sum(n[1:i])
      return(cm0[i,])
    }))
  }, uaverage = {
    cm0 <- cbind(-diag(d-1),0)
    cm <- t(sapply(1:(d-1), function(i){
      cm0[i, (i+1):d] <- n[(i+1):d]/sum(n[(i+1):d])
      return(cm0[i,])
    }))
  }, caverage = {
    cm0 <- matrix(0, nrow=d-1, ncol=d)
    cm0[,base] <- -1
    cm <- t(sapply(1:(d-1), function(i){
      id <- (1:d)[-base][1:i]
      cm0[i, id] <- n[id]/sum(n[id])
      return(cm0[i,])
    }))
  })
  colnames(cm) <- varnames
  return(cm)  
}


constrSet <- function(n, set=c("sequence","seqcontrol","lplateau","uplateau","downturn","williams"), direction=c("increase","decrease"), base=1){
  d <- length(n)
  if (d < 2) stop("less than two groups")
  if (!is.numeric(n)) stop(sQuote("n"), " is not numeric")
  if (base < 1 || base > d) stop("base is not between 1 and ", d)
  if (!is.null(names(n))) varnames <- names(n) else varnames <- 1:d
  direction <- match.arg(direction)
  switch(direction, increase = {
    direct <- 1
    vorz <- "<"
  }, decrease = {
    direct <- -1
    vorz <- ">"
  })
  set <- match.arg(set)
  switch(set, sequence = {
    sl <- lapply(1:(d-1), function(i){
      list(constr=direct*constrMat(n, type="monotone")[1:i,,drop=FALSE],
           rhs=rep(0, i),
           nec=rep(FALSE,i) )
    })
    CM <- matrix(0,nrow=1, ncol=d)
    colnames(CM) <- varnames
    sl[[d]] <- list(constr=CM, rhs=0, nec=FALSE)
    nams <- sapply(2:d, function(i) paste(varnames[1:i], collapse=vorz))
    names(sl) <- c(nams, "unconstrained")
  }, seqcontrol = {
    sl <- lapply(1:(d-1), function(i){
      list(constr=direct*constrMat(n, type="control", base=base)[1:i,,drop=FALSE],
           rhs=rep(0, i),
           nec=rep(FALSE,i) )
    })
    sl[[d]] <- list(constr=matrix(0,nrow=1, ncol=d), rhs=0, nec=FALSE)
    nams <- sapply(1:(d-1), function(i) paste(paste(varnames[base], varnames[(1:d)[-base]][1:i], sep=vorz), collapse="|"))
    names(sl) <- c(nams, "unconstrained")
  }, lplateau = {
    sl <- lapply(1:d, function(i){
      list(constr=direct*constrMat(n, type="monotone"),
           rhs=rep(0, d-1),
           nec=c(rep(TRUE,d-i), rep(FALSE,i-1)) )
    })
    sl[[d+1]] <- list(constr=matrix(0,nrow=1, ncol=d), rhs=0, nec=FALSE)
    nams <- sapply(1:d, function(i){
      vz <- rep(vorz, length=d-1)
      vz[sl[[i]]$nec] <- "="
      vz <- c(vz, "")
      paste(paste(varnames[1:d], vz, sep=""), collapse="")
    })
    names(sl) <- c(nams, "unconstrained")
  }, uplateau = {
    sl <- lapply(1:d, function(i){
      list(constr=direct*constrMat(n, type="monotone"),
           rhs=rep(0, d-1),
           nec=c(rep(FALSE,d-i), rep(TRUE,i-1)) )
    })
    sl[[d+1]] <- list(constr=matrix(0,nrow=1, ncol=d), rhs=0, nec=FALSE)
    nams <- sapply(1:d, function(i){
      vz <- rep(vorz, length=d-1)
      vz[sl[[i]]$nec] <- "="
      vz <- c(vz, "")
      paste(paste(varnames[1:d], vz, sep=""), collapse="")
    })
    names(sl) <- c(nams, "unconstrained")
  }, downturn = {
    sl <- list()
    sl[2:(d-1)] <- lapply((d-1):2, function(i){
      cm <- constrMat(n, type="monotone")
      cm[i:(d-1),] <- -1*cm[i:(d-1),]
      list(constr=cm,
           rhs=rep(0, d-1),
           nec=rep(FALSE,d-1) )
    })
    sl[[1]] <- list(constr=direct*constrMat(n, type="monotone"), rhs=rep(0, d-1),nec=rep(FALSE,d-1))
    sl[[d]] <- list(constr=matrix(0,nrow=1, ncol=d), rhs=0, nec=FALSE)
    nams <- sapply((d-1):1, function(i){
      if (vorz == ">") vz <- rep("<", length=d-1) else vz <- rep(">", length=d-1)
      vz[1:i] <- vorz
      vz <- c(vz, "")
      paste(paste(varnames[1:d], vz, sep=""), collapse="")
    })
    names(sl) <- c(nams, "unconstrained")
  }, williams = {
    sl <- lapply(1:(d-1), function(i){
      list(constr=direct*constrMat(n, type="caverage", base=1)[1:i,,drop=FALSE],
           rhs=rep(0, i),
           nec=rep(FALSE,i) )
    })
    CM <- matrix(0,nrow=1, ncol=d)
    colnames(CM) <- varnames
    sl[[d]] <- list(constr=CM, rhs=0, nec=FALSE)
    nams <- sapply(2:d, function(i) paste(varnames[1], " ", vorz, " ave(", paste(varnames[2:i], collapse=","), ")", sep=""))
    nams[1] <- paste(varnames[1], vorz, varnames[2], sep=" ")
    names(sl) <- c(nams, "unconstrained")
  })
  return(sl)
}
