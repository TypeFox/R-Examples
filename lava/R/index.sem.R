##' @export
updatelvm <- function(x,mean=TRUE,...) {
  index(x) <- reindex(x,mean=mean,...)
  x$parpos <- parpos(x,mean=mean,...)
  return(x)
}

##' @export
"index" <- function(x,...) UseMethod("index")

##' @export
"index<-" <- function(x,...,value) UseMethod("index<-")

##' @export
"index.lvm" <- function(x,...) { x$index }

##' @export
"index.lvmfit" <- function(x,...) { index(Model(x)) }

##' @export
"index<-.lvm" <- function(x,...,value)  { x$index <- value; return(x) }

##' @export
"index<-.lvmfit" <- function(x,...,value) { Model(x)$index <- value; return(x) }


###   A  ## Matrix with fixed parameters and ones where parameters are free
###   J  ## Manifest variable selection matrix
###   M0 ## Index of free regression parameters
###   M1 ## Index of free and _unique_ regression parameters
###   P  ## Matrix with fixed variance parameters and ones where parameters are free
###   P0 ## Index of free variance parameters
###   P1 ## Index of free and _unique_ regression parameters
###   npar.var  ## Number of covariance parameters
##' @export
`reindex` <-
function(x, sparse=FALSE,standard=TRUE,zeroones=FALSE,deriv=FALSE,mean=TRUE) { ## Extract indices of parameters from model
  x$parpos <- NULL
  M <- x$M

  eta <- latent(x) ## Latent variables/Factors
  m <- length(eta)
  obs <- manifest(x)  ## Manifest/Observed variables
  endo <- endogenous(x)
  exo <- exogenous(x) ##,index=FALSE)

  allvars <- vars(x)
  eta.idx <- na.omit(match(eta,allvars))
  obs.idx <- na.omit(match(obs,allvars))
  exo.idx <- na.omit(match(exo,allvars))
  exo.obsidx <- na.omit(match(exo,obs))
  endo.obsidx <- na.omit(match(endo,obs))

  fix.idx <- !is.na(x$fix) ## Index of fixed parameters
  covfix.idx <- !is.na(x$covfix) ## Index of fixed covariance parameters

  constrain.par <- NULL
  if (length(constrain(x))>0) constrain.par <- names(constrain(x))

  M0 <- M;  M0[fix.idx] <- 0 ## Matrix of indicators of free regression-parameters (removing fixed parameters)
  M1 <- M0; ## Matrix of indiciator of free _unique_ regression parameters (removing fixed _and_ duplicate parameters)
  parname <- unique(x$par[!is.na(x$par)])
##  parname.all <- unique(x$par[!is.na(x$par)])
##  parname <- setdiff(parname.all,constrain.par)
  for (p in parname) {
    ii <- which(x$par==p)
    if (length(ii)>1)
      M1[ii[-1]] <- 0
    if (p %in% constrain.par)
      M0[ii] <- M1[ii] <- 0
  }
  npar.reg <- sum(M1) ## Number of free regression parameters

  P <- x$cov;

  P0 <- P;  P0[covfix.idx] <- 0 ## Matrix of indicators of free covariance-parameters (removing fixed parameters)
  if (length(exo.idx)>0)
      P0[exo.idx,exo.idx] <- 0 ## 6/1-2011
  P1 <- P0 ## Matrix of indiciator of free _unique_ variance parameters (removing fixed _and_ duplicate parameters)
  covparname <- unique(x$covpar[!is.na(x$covpar)])
  for (p in covparname) {
    ii <- which(x$covpar==p)
    if (length(ii)>1)
      P1[ii[-1]] <- 0
    if (p%in%c(parname,constrain.par))
      P0[ii] <- P1[ii] <- 0
  }

  ##  P1. <- P1[-exo.idx,-exo.idx]
  npar.var <- sum(c(diag(P1),P1[lower.tri(P1)]))
  parnames <- paste0("p", seq_len(npar.reg+npar.var))

  A <- M
  A[fix.idx] <- x$fix[fix.idx] ## ... with fixed parameters in plac
  P[covfix.idx] <- x$covfix[covfix.idx] ## ... with fixed parameters in plac


  px <- Jy <- J <- I <- diag(nrow=length(vars(x)))
  if (m>0) {
    J[eta.idx,eta.idx] <- 0; J <- J[-eta.idx,,drop=FALSE]
  } ## Selection matrix (selecting observed variables)
  {
    ## Selection matrix (selection endogenous variables)
    if (length(c(eta.idx,exo.idx))>0) {
      Jy[c(eta.idx,exo.idx),c(eta.idx,exo.idx)] <- 0; Jy <- Jy[-c(eta.idx,exo.idx),,drop=FALSE]
    }
    ## Cancelation matrix (cancels rows with exogenous variables)
    px[exo.idx,exo.idx] <- 0
  }

  ## Creating indicitor of free mean-parameters
  fixed <- sapply(x$mean, function(y) is.numeric(y) & !is.na(y))
  named <- sapply(x$mean, function(y) is.character(y) & !is.na(y))
  mparname <- NULL
  if (length(named)>0)
      mparname <- unlist(unique(x$mean[named]))
  v0 <- rep(1,length(x$mean)) ## Vector of indicators of free mean-parameters

  v0[exo.idx] <- 0
  if (length(fixed)>0) v0[fixed] <- 0;
  v1 <- v0
  for (p in mparname) {
    idx <- which(x$mean==p)
    if (length(idx)>1) {
##      print(idx[-1])
      v1[idx[-1]] <- 0
    }
    if (p%in%c(parname,covparname,constrain.par))
      v0[idx] <- v1[idx] <- 0
  } ## duplicate parameters

  ###
  ### Extra parameters
  ###
  efixed <- sapply(x$exfix, function(y) is.numeric(y) & !is.na(y))
  enamed <- sapply(x$exfix, function(y) is.character(y) & !is.na(y))
  eparname <- unlist(unique(x$exfix[enamed]))
  ## Extra parameters
  e0 <- rep(1,length(x$expar)) ## Indicators of free extra par.
  if (length(efixed)>0)
    e0[efixed] <- 0
  e1 <- e0
  for (p in eparname) {
    idx <- which(x$exfix==p)
    if (length(idx)>1) {
      e1[idx[-1]] <- 0
    }
    if (p%in%c(parname,covparname,constrain.par,mparname))
      e0[idx] <- e1[idx] <- 0
  } ## duplicate parameters


  ## Return:
  ## Adjacency-matrix (M)
  ## Matrix of regression-parameters (0,1) _with_ fixed parameters (A)
  ## Matrix of variance-parameters (indicators 0,1) (P)
  ## Manifest selection matrix (J),
  ## Position of variables matrix (Apos),
  ## Position of covariance variables matrix (Ppos),
  ## Position/Indicator matrix of free regression parameters (M0)
  res <- list(vars=allvars, manifest=obs, exogenous=exo, latent=eta,
              endogenous=endo,
              exo.idx=exo.idx, eta.idx=eta.idx,
              exo.obsidx=exo.obsidx, endo.obsidx=endo.obsidx,
              obs.idx=obs.idx,
              endo.idx=setdiff(obs.idx,exo.idx))

  if (standard) {
    res <- c(res, list(M=M, A=A, P=P,
                       P0=P0, P1=P1,
                       M0=M0, M1=M1,
                       v0=v0, v1=v1,
                       e0=e0, e1=e1,
                       npar=(npar.reg+npar.var),
                       npar.reg=npar.reg,
                       npar.var=npar.var,
                       npar.mean=sum(v1),
                       npar.ex=sum(e1),
                       constrain.par=constrain.par))
    npar.total <- res$npar+res$npar.mean+res$npar.ex
    which.diag <- NULL
    if (length(P1)>0)
        which.diag <- which(diag(P1==1))

    res <- c(res, list(parname.all=parname, parname=setdiff(parname,constrain.par),
                       which.diag=which.diag,
                       covparname.all=covparname,
                       covparname=setdiff(covparname,constrain.par),
                       meanfixed=fixed, meannamed=named,
                       mparname.all=mparname,
                       mparname=setdiff(mparname,constrain.par),
                       eparname.all=eparname,
                       eparname=setdiff(eparname,constrain.par),
                       J=J, Jy=Jy, px=px, sparse=sparse))

    parname.all.reg.idx <- parname.all.reg.tidx <-
      parname.reg.tidx <- parname.reg.idx <- c()
    for (p in res$parname.all) {
      ipos <- which((x$par==p))
      tipos <- which(t(x$par==p))
      if (p%in%res$parname) {
        parname.reg.idx <- c(parname.reg.idx, list(ipos))
        parname.reg.tidx <- c(parname.reg.tidx, list(tipos))
      }
      parname.all.reg.idx <- c(parname.all.reg.idx, list(ipos))
      parname.all.reg.tidx <- c(parname.all.reg.tidx, list(tipos))
    };
    if (length(parname.reg.idx)>0) {
      names(parname.reg.idx) <- names(parname.reg.tidx) <- res$parname
    }
    if (length(parname.all.reg.idx)>0) {
      names(parname.all.reg.idx) <- names(parname.all.reg.tidx) <- res$parname.all
    }
    covparname.all.idx <- covparname.idx <- c()
    for (p in res$covparname.all) {
      ipos <- which(x$covpar==p)
      if (p%in%res$covparname)
        covparname.idx <- c(covparname.idx, list(ipos))
      covparname.all.idx <- c(covparname.all.idx, list(ipos))
    };
    if (length(covparname.idx)>0)
      names(covparname.idx) <- res$covparname
    if (length(covparname.all.idx)>0)
      names(covparname.all.idx) <- res$covparname.all

    mparname.all.idx <- mparname.idx <- c()
    for (p in res$mparname.all) {
      ipos <- which(x$mean==p)
      if (p%in%mparname)
        mparname.idx <- c(mparname.idx, list(ipos))
      mparname.all.idx <- c(mparname.all.idx, list(ipos))
    };
    if (length(mparname.idx)>0)
      names(mparname.idx) <- res$mparname
    if (length(mparname.all.idx)>0)
      names(mparname.all.idx) <- res$mparname.all

    eparname.all.idx <- eparname.idx <- c()
    for (p in res$eparname.all) {
      ipos <- which(x$exfix==p)
      if (p%in%eparname)
        eparname.idx <- c(eparname.idx, list(ipos))
      eparname.all.idx <- c(eparname.all.idx, list(ipos))
    };
    if (length(eparname.idx)>0)
      names(eparname.idx) <- res$eparname
    if (length(eparname.all.idx)>0)
      names(eparname.all.idx) <- res$eparname.all


    res <- c(res, list(mparname.idx=mparname.idx,
                       covparname.idx=covparname.idx,
                       parname.reg.idx=parname.reg.idx,
                       parname.reg.tidx=parname.reg.tidx,
                       mparname.all.idx=mparname.all.idx,
                       eparname.all.idx=eparname.all.idx,
                       covparname.all.idx=covparname.all.idx,
                       parname.all.reg.idx=parname.all.reg.idx,
                       parname.all.reg.tidx=parname.all.reg.tidx
                       ))

  } else {
    res <- index(x)
  }

  if (zeroones) {
    if (sparse) {
      if (!requireNamespace("Matrix",quietly=TRUE)) stop("package Matrix not available")
      Ik <- Matrix::Diagonal(length(obs))
      Im <- Matrix::Diagonal(ncol(A))
      Kkk <- NULL
      J <- as(J, "sparseMatrix")
      Jy <- as(Jy, "sparseMatrix")
      px <- as(px, "sparseMatrix")

    } else {
      Ik <- diag(nrow=length(obs))
      Im <- diag(nrow=ncol(A))
    }
    Kkk <- NULL


    res[c("Ik","Im","Kkk")] <- NULL
    res <- c(res, list(Ik=Ik, Im=Im, Kkk=Kkk))
  }
  if (deriv && length(P)>0) {
    if (res$npar.mean>0 & mean)
      D <- deriv.lvm(x,meanpar=rep(1,res$npar.mean),zeroones=TRUE)
    else
      D <- deriv.lvm(x,meanpar=NULL,zeroones=TRUE)
    res[c("dA","dP","dv")] <- NULL
    res <- c(res, list(dA=D$dA, dP=D$dP, dv=D$dv))
  }

  if (length(P)>0)
  res <- c(res,mat.lvm(x,res))

  return(res)
}
