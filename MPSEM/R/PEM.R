#
### MPSEM package R functions and C wrappers.
#
PEMInfluence <- function(x,mroot=TRUE) {
  if(attr(x,"class") != "graph")
    stop("Parameter 'x' must be of class 'graph'")
  from <- x$edge[[1L]] ; to <- x$edge[[2L]] ; ev <- attr(x,"ev")
  roots <- (1L:ev[2L])[is.na(match(1L:ev[2L],to))] ; nroots <- length(roots)
  ## If several roots, add a guide vertex and nroots edges connecting each root from the guide.
  if(nroots != 1L) {
    if(!mroot)
      stop("Multiple (",nroots,") roots are detected but multiple rooting is not allowed (i.e. mroot == FALSE).")
    warning("Multiple (",nroots,") roots detected.")
    from <- c(rep(1L,nroots),from+1L)
    to <- c(roots,to)+1L
    ev <- ev+c(nroots,1L)
  }
  B <- matrix(.C("PEMInfMat",
              as.integer(from),
              as.integer(to),
              as.integer(ev[1L]),
              as.integer(ev[2L]),
              B = integer(ev[2L]*ev[1L]))$B,ev[2L],ev[1L])
  if(nroots != 1L)
    B <- B[-1L,]    ## Strip the guidance vertex if many roots.
  if(!is.null(attr(x,"vlabel")))
    rownames(B) <- attr(x,"vlabel")
  if(!is.null(attr(x,"elabel"))) {
    if(nroots == 1L) {  ## Edge labels should be handled with 
      colnames(B) <- attr(x,"elabel")
    } else {
      colnames(B) <- c(paste("V->",attr(x,"vlabel")[roots],sep=""),attr(x,"elabel"))
    }
  }
  return(B)
}
#
PEMweights <- function (d,a=0,psi=1) {
    nd <- length(d)
    a <- rep(a, length.out = nd)
    psi <- rep(psi, length.out = nd)
    return(.C("PEMweight",
              as.double(d),
              as.integer(nd),
              as.double(a),
              as.double(psi),
              res = double(nd))$res)
}
#
PEM.build <- function(x,d="distance",sp="species",a=0,psi=1,tol=.Machine$double.eps^0.5) {
  if(attr(x,"class") != "graph")
    stop("Parameter 'x' must be of class 'graph'")
  if(is.null(x$edge[[d]]))
    stop("There is no property '",d,"' to be used as edges lengths.")
  if(is.null(x$vertex[[sp]]))
    stop("There is no property '",sp,"' to indicate species vertices.")
  nsp <- sum(x$vertex[[sp]])
  ### All graphs sould be single-rooted
  ev <- as.integer(attr(x, "ev"))  # Just to be sure.
  a <- rep(a,length.out=ev[1L])
  psi <- rep(psi,length.out=ev[1L])
  out <- list(x=x,sp=x$vertex[[sp]])
  out[["B"]] <- matrix(.C("PEMInfMat",
                        as.integer(x$edge[[1L]]),
                        as.integer(x$edge[[2L]]),
                        ev[1L],ev[2L],
                        B = integer(ev[2L]*ev[1L]))$B,ev[2L],ev[1L])[x$vertex[[sp]],]
  out <- c(out,.C("PEMbuildC",
                  ne=ev[1L],nsp=nsp,
                  Bc=as.double(out$B),
                  means=double(ev[1L]),
                  dist=as.double(x$edge[[d]]),
                  a=as.double(a),
                  psi=as.double(psi),
                  w=double(ev[1L]),
                  BcW=double(nsp*ev[1L])))
  attr(out$Bc,"dim") <- c(nsp,ev[1L])
  attr(out$BcW,"dim") <- c(nsp,ev[1L])
  dimnames(out$Bc) <- dimnames(out$BcW) <- list(attr(x,"vlabel")[x$vertex[[sp]]],attr(x,"elabel"))
  out <- c(out, La.svd(out$BcW,nsp,nsp))
  sel <- out$d >= tol
  out$d <- out$d[sel]
  out$u <- out$u[,sel,drop=FALSE]
  out$vt <- out$vt[sel,,drop=FALSE]
  rownames(out$vt) <- colnames(out$u) <- paste("V",1L:sum(sel),sep="_")
  rownames(out$u) <- attr(x,"vlabel")[x$vertex[[sp]]]
  colnames(out$vt) <- attr(x,"elabel")
  attr(out,"class") <- "PEM"
  return(out)
}
#
print.PEM <- function(x,...) {
  cat("A phylogenetic eigenvector map (PEM) for ",x$nsp," species:\n")
  if(x$nsp >= 10L)
    cat(paste(rownames(x$u)[1L:8L],collapse=","),"...",rownames(x$u)[nrow(x$u)],"\n")
  else
    cat(paste(rownames(x$u),collapse=","),"\n")
  cat("obtained from the following phylogenetic graph:\n")
  MPSEM::print.graph(x$x)
  return(invisible(NULL))
}
#
PEM.updater <- function(object,a,psi=1,tol=.Machine$double.eps^0.5) {
  nsp <- object$nsp
  ne <- object$ne
  object$a <- rep(a,length.out=ne)
  object$psi <- rep(psi,length.out=ne)
  out <- .C("PEMupdateC",
            as.integer(ne),
            as.integer(nsp),
            as.double(object$Bc),
            as.double(object$dist),
            as.double(object$a),
            as.double(object$psi),
            w=as.double(object$w),
            BcW=as.double(object$BcW))
  object$w <- out$w
  object$BcW[] <- out$BcW
  out <- La.svd(object$BcW,nsp,nsp)
  sel <- out$d > tol
  object$d <- out$d[sel]
  object$u[] <- out$u[, sel, drop = FALSE]
  object$vt[] <- out$vt[sel, , drop = FALSE]
  return(object)
}
# Coerce class PEM in data.frame. Useful, i.e., for use in lm() which calls that method on its parameter 'data'.
as.data.frame.PEM <- function(x, row.names = NULL, optional = FALSE, ...) {
  return(as.data.frame(x$u))
}
#
predict.PEM <- function (object, targets, lmobject, newdata, interval = c("none", "confidence", "prediction"), level = 0.95, ...) {
  if(missing(newdata)) newdata <- Locations2PEMscores(object, targets) else {
    if(nrow(targets$locations)!=nrow(newdata))
      stop("'newdata' has ",nrow(newdata)," rows but the number of target species is ",nrow(targets$locations),".")
    tmp <- Locations2PEMscores(object, targets)
    rownames(newdata) <- rownames(tmp$scores)
    tmp$scores <- cbind(newdata,tmp$scores)
    newdata <- tmp
    rm(tmp)
  }
  interval <- match.arg(interval)
  Residual.variance <- diag(t(lmobject$residuals) %*% lmobject$residuals)/lmobject$df
  Xh <- cbind(1, as.matrix(newdata$scores[, attr(lmobject$terms, "term.labels"), drop = FALSE]))
  pred <- Xh %*% lmobject$coefficients
  if (interval == "none") return(pred)
  R <- qr.R(lmobject$qr)
  invXtX <- solve(t(R) %*% R)
  XhinvXtXtXh <- diag(Xh %*% invXtX %*% t(Xh))
  if (interval == "confidence")
    S <- sqrt(t((newdata$VarianceFactor/nrow(object$y)) + 
         matrix(Residual.variance, length(Residual.variance), 
         length(XhinvXtXtXh)) * XhinvXtXtXh))
  if (interval == "prediction") 
    S <- sqrt(t(newdata$VarianceFactor + matrix(Residual.variance, 
         length(Residual.variance), length(XhinvXtXtXh)) * 
         (1 + XhinvXtXtXh)))
  return(list(values = pred, lower = pred + S * qt(0.5 * (1 - level), lmobject$df),
              upper = pred + S * qt(0.5 * (1 - level), lmobject$df, lower.tail = FALSE)))
}
#
PEM.fitSimple <- function(y,x,w,d="distance",sp="species",lower=0,upper=1,tol=.Machine$double.eps^0.5) {
  object <- PEM.build(w,d=d,sp=sp,a=lower)
  nsp <- sum(w$vertex[[sp]])
  if(is.matrix(y)) {
    if(nrow(y)!=nsp)
      stop("The number of rows in 'y' must match the number of species in phylogenetic graph 'w'")
    p <- ncol(y)
  } else {
    if(length(y)!=nsp)
      stop("The number of elements in 'y' must match the number of species in phylogenetic graph 'w'")
    y <- cbind(y) ; p <- 1L
  }
  f <- function(par,y,x,object,nsp,p) {
    object <- PEM.updater(object,a=par[1L],psi=1,tol=tol)
    Sinv <- object$u %*% diag(object$d^(-2)) %*% t(object$u)
    logdetS <- sum(log(object$d^2))
    if(is.null(x)) {
      BX <- mean(y)
      res <- y-BX
    } else {
      BX <- MASS::ginv(t(x)%*%Sinv%*%x) %*% t(x) %*% Sinv %*% y
      res <- y-x%*%BX
      BX <- c(mean(res),BX)
      res <- res-BX[1]
    }
    dev <- 0
    for(i in 1L:p) {
      S2 <- as.numeric(t(res[,i,drop=FALSE])%*%Sinv%*%res[,i,drop=FALSE]/nsp)
      dev <- dev + nsp + (nsp-1L)*log(2*pi) + nsp*log(S2) + logdetS
    }
    return(dev)
  }
  opt <- optim(par=0,f,method="L-BFGS-B",lower=lower,upper=upper,y=y,x=x,object=object,nsp=nsp,p=p)
  if(opt$convergence)
    warning("No optimum found... Message from optim() - ",opt$message,". Status = ",opt$convergence)
  object <- PEM.updater(object,a=opt$par[1L],psi=1,tol=tol)
  Sinv <- object$u %*% diag(object$d^(-2)) %*% t(object$u)
  if(is.null(x)) {
    BX <- mean(y)
    res <- y-BX
  } else {
    BX <- MASS::ginv(t(x)%*%Sinv%*%x) %*% t(x) %*% Sinv %*% y
    res <- y-x%*%BX
    BX <- c(mean(res),BX)
    res <- res-BX[1]
  }
  S2 <- numeric(p)
  for(i in 1L:p)
    S2[i] <- as.numeric(t(res[,i,drop=FALSE])%*%Sinv%*%res[,i,drop=FALSE]/nsp)
  object$S2 <- S2
  object$y <- y
  object$optim <- opt
  object$x$edge[["a"]] <- rep(opt$par[1L],attr(object$x,"ev")[1L])
  return(object)
}
#
PEM.forcedSimple <- function(y,x,w,d="distance",sp="species",a=0,psi=1,tol=.Machine$double.eps^0.5) {
  object <- PEM.build(w,d=d,sp=sp,a=a,psi=psi,tol=tol)
  nsp <- sum(w$vertex[[sp]])
  if(is.matrix(y)) {
    if(nrow(y)!=nsp)
      stop("The number of rows in 'y' must match the number of species in phylogenetic graph 'w'")
    p <- ncol(y)
  } else {
    if(length(y)!=nsp)
      stop("The number of elements in 'y' must match the number of species in phylogenetic graph 'w'")
    y <- cbind(y) ; p <- 1L
  }
  Sinv <- object$u %*% diag(object$d^(-2)) %*% t(object$u)
  if(is.null(x)) {
    BX <- mean(y) ; res <- y-BX
  } else {
    BX <- MASS::ginv(t(x)%*%Sinv%*%x) %*% t(x) %*% Sinv %*% y
    res <- y-x%*%BX
    BX <- c(mean(res),BX)
    res <- res-BX[1]
  }
  S2 <- numeric(p)
  for(i in 1L:p)
    S2[i] <- as.numeric(t(res[,i,drop=FALSE])%*%Sinv%*%res[,i,drop=FALSE]/nsp)
  object$S2 <- S2
  object$y <- y
  object$optim <- list()
  object$optim$par <- c(a,psi)
  object$x$edge[["a"]] <- rep(a,attr(object$x,"ev")[1L])  # Hence the 'Simple' in 'fitSimple'
  object$x$edge[["psi"]] <- rep(psi,attr(object$x,"ev")[1L])
  return(object)
}
#
getGraphLocations <- function(tpall,targets) {
  oldnlab <- tpall$node.label
  tpall$node.label <- newnlab <- paste("N",1:tpall$Nnode,sep="")
  tpmodel <- drop.tip(tpall, targets) ; tpmodel$root.edge <- tpall$root.edge
  xmodel <- Phylo2DirectedGraph(tpmodel)
  Bmodel <- PEMInfluence(xmodel)
  loc <- matrix(NA, length(targets), ncol(Bmodel), dimnames = list(targets,colnames(Bmodel)))
  dtt <- rep(NA, length(targets))
  for (i in 1L:length(targets)) {    # Target species need to be striped one-by-one.
    tptargetonly <- if (length(targets) == 1L) tpall else drop.tip(tpall, targets[-i]) # plot(tptargetonly,show.node.label=TRUE)
    tptargetonly$root.edge <- tpall$root.edge
    Vnames <- c(tptargetonly$tip.label,tptargetonly$node.label)              # The vertices names.
    VBidx <- which(tptargetonly$edge[,2L] == which(Vnames == targets[i]))    # The index of the vertex immediatly below.
    VB <- tptargetonly$edge[VBidx,1L]                                        # VB: the vertex immediatly below
    VBName <- Vnames[VB]                                                     # The name of VB
    dBX <- tptargetonly$edge.length[VBidx]                                   # Distance between VB and VX
    VA <- tptargetonly$edge[tptargetonly$edge[,2L] == VB,1L]                 # VA: the vertex below the one immediatly below and
    VAName <- Vnames[VA]                                                            # its name. That vertex may not exist.
    dAB <- tptargetonly$edge.length[tptargetonly$edge[,2L] == VB]                   # Distance between VB and VA (if exists)
    VC <- tptargetonly$edge[tptargetonly$edge[,1L] == VB,2L]                        # The vertices above the vertex immediatly below and
    VCName <- Vnames[VC]                                                            # their names. 2 or more (in case of tri+ chotomy).
    dBC <- tptargetonly$edge.length[tptargetonly$edge[,1L] == VB]                         # Distances above VB.
    VC <- VC[VCName!=targets[i]] ; dBC <- dBC[VCName!=targets[i]] ; VCName <- Vnames[VC]  # Strip the target itself from the VC list.
    if(length(VA)==0) {                                                                   # If the target species is beyond the root:
      loc[i,] <- 0                                                      # In all case: location = the root.
      dtt[i] <- dBX+min(dBC)                                            # Distance between the off-root species and the remaining ones.
    } else {                                                            # If the target species is not beyond the root:
      dtt[i] <- dBX                                          # In all cases: dtt == the distance between the target and the vertex immediatly below.
      if(length(VC)>1) {                                     # When VB is more than dichotomic (several Vertices C):
        loc[i,] <- Bmodel[VBName,] * xmodel$edge$distance    # Coordinates are those of the still-existing vertex B.
      } else {                                               # When VB collapses (it was dichotomic):
        loc[i,] <- Bmodel[VAName,] * xmodel$edge$distance    # Copy the coordinates of vertex VA.
        loc[i,Bmodel[VAName,] != Bmodel[VCName,]] <- dAB     # Change the length of that edge connecting VA to the single VC (0 for VA) to dAB.
      }
    }
  }
  if(!is.null(oldnlab)) {
    oldnlab <- oldnlab[match(attr(xmodel,"vlabel"),newnlab)]
    attr(xmodel,"vlabel")[!is.na(oldnlab)] <- na.omit(oldnlab)
  }
  return(list(x = xmodel, locations = loc, LCA2target = dtt))
}
#
getAncGraphLocations <- function(x, tpall) {
  if(missing(x) && !missing(tpall))
    x <- Phylo2DirectedGraph(tpall)
  loc <- PEMInfluence(x)[!x$vertex$species,] * x$edge$distance
  return(list(x = x, locations = loc, LCA2target = numeric(nrow(loc))))
}
#
Locations2PEMscores <- function(object, gsc) {
  if(object$ne != ncol(gsc$locations))
    stop("Numbers of edge coordinates mismatch: gsc$locations has ",ncol(gsc$locations)," columns while the phylogenetic graph has ",object$ne," edges.")
  if(is.null(object$y)||is.null(object$optim))
    stop("PEM has no attached data/optimization parameters.")
  ntgt <- nrow(gsc$locations) ; nsv <- length(object$d)
  out <- list()
  out[["scores"]] <- matrix(.C("PEMLoc2Scores",
                               as.integer(object$ne),
                               as.double(object$means*object$w),
                               as.integer(ntgt),
                               as.double(gsc$locations),
                               as.double(object$a),
                               as.double(object$psi),
                               as.integer(nsv),
                               object$d,
                               object$vt,
                               scores=double(nsv*ntgt))$scores,ntgt,nsv)
  dimnames(out[["scores"]]) <- list(rownames(gsc$locations),colnames(object$u))
  VarianceFactor <- MPSEM::PEMvar(gsc$LCA2target,a=object$optim$par[1L],psi=1)
  VarianceFactor <- matrix(object$S2,length(object$S2),1L)%*%matrix(VarianceFactor,1L,length(VarianceFactor))
  dimnames(VarianceFactor) <- list(colnames(object$y),rownames(gsc$locations))
  out[["VarianceFactor"]] <- VarianceFactor
  return(out)
}
#
