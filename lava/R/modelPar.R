
##' @export
`modelPar` <-
  function(x,p,...) UseMethod("modelPar")

###{{{ modelPar.lvmfit
##' @export
modelPar.lvmfit <- function(x, p=pars(x), ...) modelPar(Model(x),p=p,...)

###}}} modelPar.lvmfit

###{{{ modelPar

##' @export
modelPar.lvm <- function(x,p, ...) {
  npar <- index(x)$npar
  npar.mean <- index(x)$npar.mean
  if (length(p)!=npar & length(p)<(npar+npar.mean)) stop("Wrong dimension of parameter vector!")
  p2 <- NULL
  if (length(p)!=npar) { ## if meanstructure
    meanpar <- p[seq_len(npar.mean)]
    p. <- p
    if (length(meanpar)>0) {
        p. <- p[-seq_len(npar.mean)]
    } else meanpar <- NULL
    p <- p.[seq_len(npar)]
    if (npar>0) {
        p2 <- p.[-seq_len(npar)]
    } else p2 <- p.
  } else {
    meanpar <- NULL
    p2 <- NULL
  }
  return(list(p=p,meanpar=meanpar,p2=p2))
}

###}}} modelpar.lvm

###{{{ modelPar.multigroupfit

##' @export
modelPar.multigroupfit <- function(x,p=pars(x),...) {
  modelPar(Model(x),p,...)
}
###}}}

###{{{ modelPar.multigroup

##' @export
modelPar.multigroup <- function(x,p, ...) {
  if (length(p)==x$npar) {
    pp <- lapply(x$parposN,function(z) p[z])
    res <- list(p=pp, par=pp, mean=NULL)
    return(res)
  }
  Nmean <- unlist(lapply(x$meanposN,length))
  Npar <- unlist(lapply(x$parposN,length))
  ##ppos <- mapply("+",x$parposN,as.list(Nmean),SIMPLIFY=FALSE)
  ppos <- x$parposN
  pp <- lapply(ppos,function(z) p[z+x$npar.mean])

  if (length(pp)==0) pp <- lapply(seq_len(x$ngroup),function(x) logical())
  mm <- lapply(x$meanposN,function(x) p[x])
  if (is.null(mm)) mm <- lapply(seq_len(x$ngroup),logical())
  pm <- mm
  for (i in seq_len(length(pm))) pm[[i]] <- c(pm[[i]],pp[[i]])
  res <- list(p=pm,par=pp,mean=mm)
  return(res)
}

###}}}

modelPar2.multigroup <-
  function(x,p, ...) {
  npar <- x$npar
  npar.mean <- x$npar.mean
  k <- x$ngroup
  if (length(p)!=npar & length(p)!=(npar+npar.mean)) stop("Wrong dimension of parameter vector!")
  if (length(p)!=npar) { ## if meanstructure
      meanpar <- p[seq_len(npar.mean)]
      p. <- p[-seq_len(npar.mean)]
    } else {
      meanpar <- NULL
      p. <- p
    }


  parlist <- list(); for (i in seq_len(k)) parlist[[i]] <- numeric(length(x$parlist[[i]]))
  if (!is.null(meanpar)) {
    meanlist <- list(); for (i in seq_len(k)) meanlist[[i]] <- numeric(length(x$meanlist[[i]]))
  }

  if (length(p.)>0)
  for (i in seq_along(p.)) {
    for (j in seq_len(k)) {
      idx <- match(paste0("p",i), x$parlist[[j]])
      if (!is.na(idx))
        parlist[[j]][idx] <- p.[i]
      if (!is.null(meanpar)) {
        midx <- match(paste0("p",i), x$meanlist[[j]])
        if (!is.na(midx))
          meanlist[[j]][midx] <- p.[i]
      }
    }
  }

  if (!is.null(meanpar)) {
    for (i in seq_along(meanpar)) {
      for (j in seq_len(k)) {
        idx <- match(paste0("m",i), x$meanlist[[j]])
        if (!is.na(idx))
          meanlist[[j]][idx] <- meanpar[i]
      }
    }
  } else {
    meanlist <- NULL
  }
  p0 <- parlist
  for (i in seq_along(p0))
    p0[[i]] <- c(meanlist[[i]],parlist[[i]])
  return(list(p=p0, par=parlist, mean=meanlist))
}
