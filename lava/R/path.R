##' Extract all possible paths from one variable to another connected component
##' in a latent variable model. In an estimated model the effect size is
##' decomposed into direct, indirect and total effects including approximate
##' standard errors.
##'
##' @title Extract pathways in model graph
##' @export
##' @aliases path effects path.lvm effects.lvmfit
##' totaleffects
##' @seealso \code{children}, \code{parents}
##' @return If \code{object} is of class \code{lvmfit} a list with the following
##' elements is returned \item{idx}{ A list where each element defines a
##' possible pathway via a integer vector indicating the index of the visited
##' nodes. } \item{V }{ A List of covariance matrices for each path. }
##' \item{coef }{A list of parameters estimates for each path} \item{path }{A
##' list where each element defines a possible pathway via a character vector
##' naming the visited nodes in order.  } \item{edges }{Description of 'comp2'}
##'
##' If \code{object} is of class \code{lvm} only the \code{path} element will be
##' returned.
##'
##' The \code{effects} method returns an object of class \code{effects}.
##' @note For a \code{lvmfit}-object the parameters estimates and their
##' corresponding covariance matrix are also returned.  The
##' \code{effects}-function additionally calculates the total and indirect
##' effects with approximate standard errors
##' @author Klaus K. Holst
##' @keywords methods models graphs
##' @examples
##'
##' m <- lvm(c(y1,y2,y3)~eta)
##' regression(m) <- y2~x1
##' latent(m) <- ~eta
##' regression(m) <- eta~x1+x2
##' d <- sim(m,500)
##' e <- estimate(m,d)
##'
##' path(Model(e),y2~x1)
##' parents(Model(e), ~y2)
##' children(Model(e), ~x2)
##' children(Model(e), ~x2+eta)
##' effects(e,y2~x1)
##'
##' @usage
##' \method{path}{lvm} (object, to = NULL, from, ...)
##' \method{effects}{lvmfit} (object, to, from, silent=FALSE, ...)
##' @param object Model object (\code{lvm})
##' @param to Outcome variable (string). Alternatively a formula specifying
##' response and predictor in which case the argument \code{from} is ignored.
##' @param from Response variable (string), not necessarily directly affected by
##' \code{to}.
##' @param silent Logical variable which indicates whether messages are turned
##' on/off.
##' @param \dots Additional arguments to be passed to the low level functions
##' @export
path <- function(object,...) UseMethod("path")

##' @export
path.lvmfit <- function(object,to=NULL,from,...) {
  mypath <- pathM(Model(object)$M,to,from,...)
  cc <- coef(object,level=9,labels=FALSE) ## All parameters (fixed and variable)

  #cc0 <- coef(object,level=1) ## Estimated parameters
  cc0 <- coef(object,level=2) ## Estimated parameters
  i1 <- na.omit(match(rownames(cc),rownames(cc0)))
  idx.cc0 <-  which(rownames(cc)%in%rownames(cc0)); ## Position of estimated parameters among all parameters
  S <- matrix(0,nrow(cc),nrow(cc)); rownames(S) <- colnames(S) <- rownames(cc)
  V <- object$vcov
  npar.mean <- index(object)$npar.mean
#  if (object$control$meanstructure & npar.mean>0)
#    V <- V[-c(seq_len(npar.mean)),-c(seq_len(npar.mean))]
  S[idx.cc0,idx.cc0] <- V[i1,i1]  ## "Covariance matrix" of all parameters

  idx <- list()
  coefs <- list()
  V <- list()
  for (i in seq_along(mypath)) {
    xx <- mypath[[i]]
    ii <- c()
    for (j in seq_len(length(xx)-1)) {
      st <- paste0(xx[j+1], lava.options()$symbol[1], xx[j])
      ii <- c(ii, match(st,rownames(cc)))
    }
    idx <- c(idx, list(ii))
    V <- c(V, list(S[ii,ii]))
    coefs <- c(coefs, list(cc[ii]))
  }

  edges <- list()
  for (i in seq_along(mypath)) {
    p0 <- mypath[[i]]
    ee <- c()
    for (i in seq_len(length(p0)-1)) {
      ee <- c(ee, paste(p0[i],p0[i+1],sep="~"))
    }
    edges <- c(edges, list(ee))
  }
  res <- list(idx=idx,V=V,coef=coefs, path=mypath, edges=edges)
  return(res)
}

##' @export
path.lvm <- function(object,to=NULL,from,...) pathM(object$M,to=to,from=from,...)

##' @export
path.graphNEL <- function(object,to,from,...) {
  if (inherits(to,"formula")) {
    fvar <- extractvar(to)
    if (length(fvar$x)==1 & length(fvar$y)==1)
      return(path(object,to=fvar$y,from=fvar$x))
    res <- list()
    for (y in fvar$y) {
      for (x in fvar$x) {
        cat("x=",x, " y=",y, "\n")
        res <- c(res, list(path(object,to=y,from=x)))
      }
    }
    return(res)
  }
  ff <- function(g,from=1,to=NULL,res=list()) {
    M <- graph::edgeMatrix(g)
    i1 <- which(M[1,]==from)
    for (i in i1) {
      e <- M[,i]; newto <- e[2];
      if (is.null(to) || M[2,i]==to) {
        res <- c(res, list(M[,i]))
      }
      newpath <- ff(g,from=newto,to=to,list())
      if (length(newpath)>0)
      for (j in seq_along(newpath)) {
        if (is.null(to) || (tail(newpath[[j]],1)==to))
          res <- c(res, list(c(M[,i],newpath[[j]][-1])))
      }
    }
    return(res)
  }
  idxfrom <- ifelse(is.numeric(from),from,which(from==graph::nodes(object)))
  ##M <- as(object,"matrix")
  ##reachable <- acc(M,graph::nodes(object)[idxfrom])
  reachable <- graph::acc(object,graph::nodes(object)[idxfrom])[[1]]
  
  if (is.null(to)) {
    idxto <- reachable
  } else {
    idxto <- ifelse(is.numeric(to),to,which(to==graph::nodes(object)))
  }

  if (!(graph::nodes(object)[idxto] %in% names(reachable)))
##    return(structure(NULL,to=to[1],from=from[1]))
    return(NULL)
  ##    stop("No directional relationship between variables")

  mypaths <- ff(object,idxfrom,idxto)
  res <- list()
  for (i in seq_along(mypaths)) {
    res <- c(res, list(graph::nodes(object)[mypaths[[i]]]))
  }
  return(res)
}


pathM <- function(M,to,from,...) {
  nn <- colnames(M)
  if (inherits(to,"formula")) {
    fvar <- extractvar(to)
    if (length(fvar$x)==1 & length(fvar$y)==1)
      return(pathM(M,to=fvar$y,from=fvar$x))
    res <- list()
    for (y in fvar$y) {
      for (x in fvar$x) {
        cat("x=",x, " y=",y, "\n")
        res <- c(res, list(pathM(M,to=y,from=x)))
      }
    }
    return(res)
  }

  ff <- function(g,from=1,to=NULL,res=list()) {
    i1 <- which(M[from,]==1)
    for (i in i1) {
      ##      e <- M[,i]; newto <- e[2];
      if (is.null(to) || i==to) {
        res <- c(res, list(c(from,i)))
      }
      newpath <- ff(g,from=i,to=to,list())
      if (length(newpath)>0)
      for (j in seq_along(newpath)) {
        if (is.null(to) || (tail(newpath[[j]],1)==to))
          res <- c(res, list(c(c(from,i),newpath[[j]][-1])))
      }
    }
    return(res)
  }
  idxfrom <- ifelse(is.numeric(from),from,which(from==nn))
  reachable <- acc(M,nn[idxfrom])
  
  if (is.null(to)) {
    idxto <- reachable
  } else {
    idxto <- ifelse(is.numeric(to),to,which(to==nn))
  }

  if (!(nn[idxto] %in% reachable))
    return(NULL)
  ##    stop("No directional relationship between variables")

  mypaths <- ff(M,idxfrom,idxto)
  res <- list()
  for (i in seq_along(mypaths)) {
    res <- c(res, list(nn[mypaths[[i]]]))
  }
  return(res)
}
