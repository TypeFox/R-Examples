
##' Generic method for finding indeces of model parameters
##'
##' @title Generic method for finding indeces of model parameters
##' @param x Model object
##' @param \dots Additional arguments
##' @author Klaus K. Holst
##' @export
`parpos` <-
  function(x,...) UseMethod("parpos")

##' @export
parpos.default <- function(x,p,...) {
  if (is.numeric(p)) return(p)
  na.omit(match(coef(x),p))
}

##' @export
parpos.multigroup <- function(x,p,mean=TRUE,...) {
  if (missing(p)) {
    p <- unique(unlist(lapply(x$lvm, function(z) setdiff(parlabels(z),names(constrain(z))) )))
  }
  if (!is.character(p)) p <- names(p)
  p0 <- rep(NA,with(x,npar+npar.mean));
  names(p0) <- c(x$mean,x$par)
  for (i in seq_along(x$lvm)) {
    cur <- parpos(x$lvm[[i]],p=p)
    if (length(cur)>0) {
      p0[c(x$meanpos[[i]],x$parpos[[i]])[cur]] <- names(cur)
      M <- na.omit(match(names(cur),p))
      if (length(M)>0)
        p <- p[-M]
    }
    if (length(p)==0) break;
  }
  p1 <- which(!is.na(match(x$name,p)))
  p0[p1] <- x$name[p1]
  return(structure(which(!is.na(p0)),name=p0))
##  return(p0)
}

##' @export
parpos.multigroupfit <- function(x,...) parpos.multigroup(x$model0,...)

##' @export
parpos.lvm <- function(x,p,mean=TRUE,...) {
  if (!missing(p)) {
    if (!is.character(p)) p <- names(p)
    cc1 <- coef(Model(x),mean=mean,fix=FALSE)
    cc2 <- coef(Model(x),mean=mean,fix=FALSE,labels=TRUE)
    idx1 <- na.omit(match(p,cc1))
    idx11 <- na.omit(match(p,cc2))
    res <- (union(idx1,idx11));
    if (length(res)!=length(p)) {
      names(res) <- cc1[res]
    } else {
      names(res) <- p
    }
    ##    res <- idx1; res[!is.na(idx11)] <- idx11[!is.na(idx11)]
    ##    names(res) <- p
    ord <- order(res)
    res <- sort(res)
    attributes(res)$ord <- ord
    return(res)
  }
  if (mean)
    nn <- with(index(x),matrices2(x,seq_len(npar+npar.mean+npar.ex))) ## Position of parameters
  else nn <- with(index(x),matrices(x,seq_len(npar),NULL,seq_len(npar.ex)+npar))
  nn$A[index(x)$M0!=1] <- 0
  nn$P[index(x)$P0!=1] <- 0
  nn$v[index(x)$v0!=1] <- 0
  nn$e[index(x)$e0!=1] <- 0
  nn
}

##' @export
parpos.lvmfit <- parpos.lvm
