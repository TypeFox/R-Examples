
##' @export
totaleffects <- function(object,...,value) UseMethod("totaleffects")

##' @export
totaleffects.lvmfit <- function(object,to,...,level=0.95) {
  p <- (1-level)/2
  q <- qnorm(p)
  res <- c()
  if (inherits(to,"formula")) {
    if (substr(deparse(to[3]),1,1)==".") {
      trim <- function(x) sapply(x,function(z) gsub(" ","",z,fixed=TRUE))
      to <- trim(strsplit(deparse(to),"~")[[1]][1])
    } else {
      to <- list(to)
    }
  }
  if (is.null(list(...)$from) & is.character(to)[1]) {
    to <- lapply(paste(to,setdiff(vars(object),to),sep="~"),FUN=as.formula)
  }
  ef <- function(tt) {
    f <- effects(object,tt,...)
    rbind(with(f$totalef,c(est,sd,est/sd,2*(pnorm(abs(est/sd),lower.tail=FALSE)),est+q*sd,est-q*sd)))
  }
  if (is.list(to)) {
    for (tt in to) {
      res <- rbind(res,ef(tt))
    }
  }
  else
    res <- ef(to)
  colnames(res) <- c("Estimate","Std.Err","z value","Pr(>|z|)",
                     paste0(c(1-p,p)*100,"%"))
  rownames(res) <- to
  res
}

##' @export
effects.lvmfit <- function(object,to,from,silent=FALSE,...) {
  if (missing(to)) {
    return(summary(object))
  }
  P <- path(object,to=to,from=from,...)
  if (is.null(P$path)) {
    if (inherits(to,"formula")) {
      f <- extractvar(to)
      to <- f$y; from <- f$x
    }
  } else {
    from <- P$path[[1]][1]
    to <- tail(P$path[[1]],1)
  }
  cc <- coef(object,level=9,labels=FALSE) ## All parameters (fixed and variable)
  cc0 <- cbind(coef(object)) ## Estimated parameters
  i1 <- na.omit(match(rownames(cc),rownames(cc0)))
  idx.cc0 <-  which(rownames(cc)%in%rownames(cc0)); ## Position of estimated parameters among all parameters
  S <- matrix(0,nrow(cc),nrow(cc)); rownames(S) <- colnames(S) <- rownames(cc)
  V <- object$vcov
  npar.mean <- index(object)$npar.mean
##  if (object$control$meanstructure & npar.mean>0)
##    V <- V[-seq_len(npar.mean),-seq_len(npar.mean)]
  S[idx.cc0,idx.cc0] <- V[i1,i1] ## "Covariance matrix" of all parameters

  cclab <- rownames(coef(object,level=9,labels=TRUE)) ## Identify equivalence constraints
  cctab <- table(cclab)
  equiv <- which(cctab>1)
  for (i in seq_len(length(equiv))) {
    orgpos <- which(cclab==(names(equiv)[i]))
    pos <- orgpos[-1]
    for (p in pos)
        S[p,-orgpos[1]] <- S[-orgpos[1],p] <- S[orgpos[1],-p]
  }

  idx.orig <- unique(unlist(P$idx))
  coefs.all <- cc[idx.orig]

  S.all <- S[idx.orig,idx.orig]
  idx.all <- numberdup(unlist(P$idx))
  pos <- 1; idx.list <- P$idx; for (i in seq_len(length(idx.list))) {
    K <- length(idx.list[[i]])
    idx.list[[i]] <- idx.all[pos:(pos+K-1)]; pos <- pos+K
  }
  margef <- list()
  if (length(coefs.all)==1 && is.na(coefs.all)) {
    totalef <- list(est=0,sd=0)
    margef <- c(margef,list(est=0,sd=NA))
  } else {
    totalef <- prodsumdelta(coefs.all, idx.list, S.all,...)
    for (i in seq_len(length(idx.list))) {
      margef <- c(margef, list(prodsumdelta(coefs.all, idx.list[i], S.all,...)))
    }
    paths <- list()
  }

  directidx <- which(lapply(P$path,length)==2)

  inef.list <- idx.list
  if (length(directidx)==0) {
    directef <- list(est=0, sd=NA)
  } else {
    inef.list <- inef.list[-directidx]
    directef <- margef[[directidx]]
  }
  if (length(inef.list)==0) {
    totalinef <- list(est=0,sd=NA,grad=NA,hess=NA)
  } else {
    totalinef <- prodsumdelta(coefs.all, inef.list, S.all,...)
  }

  val <- list(paths=P$path, totalef=totalef, directef=directef, totalinef=totalinef, margef=margef, from=from, to=to)
  class(val) <- "effects"
  ##    res <- c(res, list(val))
  val
}

##' @export
print.effects <- function(x,...) {
  with(x, {
    cat("\nTotal effect of '", from, "' on '", to, "':\n", sep="")
    cat("\t\t", totalef$est, " (Approx. Std.Err = ", totalef$sd, ")\n", sep="")
    cat("Direct effect of '", from, "' on '", to, "':\n", sep="")
    cat("\t\t", directef$est, " (Approx. Std.Err = ", directef$sd, ")\n", sep="")
    cat("Total indirect effect of '", from, "' on '", to, "':\n", sep="")
    cat("\t\t", totalinef$est, " (Approx. Std.Err = ", totalinef$sd, ")\n", sep="")

    cat("Indirect effects:\n");
    for (i in seq_len(length(margef))) {
      if (length(paths[[i]])>2) {
        cat("\tEffect of '", from, "' via ", paste(paths[[i]],collapse="->"), ":\n", sep="");
        cat("\t\t", margef[[i]]$est, " (Approx. Std.Err = ", margef[[i]]$sd, ")\n", sep="")
      }
    }
  })
  cat("\n");
  invisible(x)
}

##' @export
coef.effects <- function(object,...) {
  totalef <- with(object$totalef, cbind(est,sd[1]))
  directef <- with(object$directef, cbind(est,sd[1]))
  totindirectef <- with(object$totalinef, cbind(est,sd[1]))
  rownames(totalef) <- "Total"
  rownames(directef) <- "Direct"
  rownames(totindirectef) <- "Indirect"
  nn <- indirectef <- c()
  K <- seq_len(length(object$margef))
  for (i in K) {
    if (length(object$paths[[i]])>2) {
      nn <- c(nn,paste(rev(object$paths[[i]]),collapse=lava.options()$symbol[1]))
      indirectef <- rbind(indirectef, with(object$margef[[i]], c(est,sd)))
      }
  }; rownames(indirectef) <- nn
  mycoef <- rbind(totalef,directef,totindirectef,indirectef)
  mycoef <- cbind(mycoef,mycoef[,1]/mycoef[,2])
  mycoef <- cbind(mycoef,2*(pnorm(abs(mycoef[,3]),lower.tail=FALSE)))
  colnames(mycoef) <- c("Estimate","Std.Err","z value","Pr(>|z|)")
  mycoef
}

##' @export
confint.effects <- function(object,parm,level=0.95,...) {
  mycoef <- coef(object)
  p <- 1-(1-level)/2
  res <- mycoef[,1] +  + qnorm(p)*cbind(-1,1)%x%mycoef[,2]
  colnames(res) <- paste0(c(1-p,p)*100,"%")
  rownames(res) <- rownames(mycoef)
  res
}


prodtrans <- function(betas) {
  k <- length(betas)
  res <- prod(betas)
  ##  if (all(betas>0)) {
  ##    attr(res,"gradient") <- res/betas
  ##    return(res)
  ##  }
  nabla <- numeric(k)
  for (i in seq_len(k))
    nabla[i] <- prod(betas[-i])

  H <- matrix(0,k,k)
  if (k>1)
    for (i in seq_len(k-1))
      for (j in (i+1):k)
        H[j,i] <- H[i,j] <- prod(c(1,betas[-c(i,j)]))
  attr(res,"gradient") <- nabla
  attr(res,"hessian") <- H
  return(res)
}
prodsumdelta <- function(betas,prodidx,S,order=1) { ## Delta-method
  k <- length(prodidx)
  p <- length(betas)
  if (p==1) {
    return(list(est=betas, sd=sqrt(S), grad=0, hess=0))
  }
  val <- 0; grad <- numeric(p)
  H <- matrix(0,p,p)
  for (i in seq_len(k)) {
    ii <- prodidx[[i]]
    myterm <- prodtrans(betas[ii]);
    if (order>1) {
      H0 <- attributes(myterm)$hessian
      Sigma <- S[ii,ii]
       ## print(Sigma)
       ## print(H0)
       ## print(Sigma%*%H0)
      print(sum(diag(Sigma%*%H0))/2)
      val <- val + (myterm + sum(diag(Sigma%*%H0))/2)
    } else {
      val <- val + myterm
    }
    grad[ii] <- grad[ii] + attributes(myterm)$gradient
  }; grad <- matrix(grad,ncol=1)
  return(list(est=val, sd=sqrt(t(grad)%*%S%*%grad), grad=grad, hess=H))
}
