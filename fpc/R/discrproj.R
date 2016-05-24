tdecomp <- function(m){
  wm <- eigen(m, symmetric=TRUE)
  p <- ncol(m)
  wmd <- wm$values
  for (i in 1:p){
    if (abs(wmd[i])<1e-6)
      wmd[i] <- 1e-6
  }
  out <- t(wm$vectors %*% diag(sqrt(wmd)))
  out
}

discrcoord <- function(xd, clvecd, pool="n", ...) {
  x <- as.matrix(xd)
  clvec <- as.integer(clvecd)
  n <- nrow(x)
  p <- ncol(x)
  clf <- factor(clvec)
  cll <- as.integer(levels(clf)) 
  clnum <- length(cll)
  cln <- rep(0, times=clnum)
  for (i in 1:clnum){
    cln[i] <- sum(clvec==cll[i])
  }
  W <- rep(0, times=p*p)
  dim(W) <- c(p,p)
  for (i in 1:clnum){
    clx <- rep(0, times=p*cln[i])
    dim(clx) <- c(cln[i],p)
    for (j in 1:p){
      clx[,j] <- x[,j][clvec==cll[i]]
    }
    cclx <- cov(clx)
    if (cln[i]<2) cclx <- 0
    if (pool=="n")
      W <- W + ((cln[i]-1)*cclx)
    else
      W <- W + (n-1)*cclx/clnum
  }
  Tm <- tdecomp(W)
  Tinv <- solve(Tm)
  S <- (n-1)*cov(x)
  B <- S-W
  Z <- t(Tinv) %*% B %*% Tinv
  dc <- eigen(Z, symmetric=TRUE)
  units <- Tinv %*% dc$vectors * sqrt(n-clnum)
  proj <- x %*% units    
  list(ev=dc$values, units=units, proj=proj, W=W)
}

batvarcoord <- function(xd, clvecd, clnum=1){
  x <- as.matrix(xd)
  clvec <- as.integer(as.integer(clvecd)==as.integer(clnum))
  n <- nrow(x)
  p <- ncol(x)
  cll <- c(0,1)
  cln <- rep(0, times=2)
  for (i in 1:2)
    cln[i] <- sum(clvec==cll[i])
  clx <- list()
  for (i in 1:2){
    clx[[i]] <- rep(0, times=p*cln[i])
    dim(clx[[i]]) <- c(cln[i],p)
    for (j in 1:p)
      clx[[i]][,j] <- x[,j][clvec==cll[i]]
  }
  S1 <- cov(clx[[1]])
  S2 <- cov(clx[[2]])
  W <- solve(S2) %*% S1
  Weigen <- eigen(W)
  rev <- Weigen$values + 1/Weigen$values + 2
  dw <- diag(t(Weigen$vectors) %*% S2 %*% Weigen$vectors)
  svw <- matrix(ncol=p, nrow=p)
  for (i in 1:p)
    svw[,i] <- Weigen$vectors[,i]/sqrt(dw[i])
  units <- svw[,(order(-rev))]
  proj <- x %*% units    
  list(ev=Weigen$values[order(-rev)], rev=rev[(order(-rev))], 
       units=units, proj=proj, W=W, S1=S1, S2=S2)
}

batcoord <- function(xd, clvecd, clnum=1, dom="mean"){
  x <- as.matrix(xd)
  clvec <- as.integer(clvecd)
  clf <- factor(clvec)
  cll <- as.integer(levels(clf)) 
  clz <- length(cll)
  p <- ncol(x)
  if (clz!=2){
    clvec <- as.integer(clvecd==clnum)
    cll <- c(0,1)
    print("Cluster indicator has more than 2 values")
  }
  if (dom=="mean"){
    dcx <- discrcoord(x, clvec, pool="equal")
    x2 <- dcx$proj[,2:p]
    batx <- batvarcoord(x2, as.integer(clvecd==clnum))
    ev <- c(dcx$ev[1],batx$ev)
    rev <- c(max(batx$rev)+1, batx$rev)
    units <- cbind(dcx$units[,1],dcx$units[,2:p] %*% batx$units)
    proj <- cbind(dcx$proj[,1],batx$proj)
  }
  else{
    batx <- batvarcoord(x, as.integer(clvecd==clnum))
    cln <- rep(0, times=2)
    mx <- matrix(nrow=nrow(x),ncol=p)
    for (i in 1:2)
      cln[i] <- sum(clvec==cll[i])
    for (i in 1:2){
      clx <- rep(0, times=p*cln[i])
      dim(clx) <- c(cln[i],p)
      for (j in 1:p){
        clx[,j] <- x[,j][clvec==cll[i]]
      }
      mx[i,] <- colMeans(clx)
    }
    mdiff <- mx[2,]-mx[1,]
    ev <- batx$ev
    rev <- rep(0, times=p)
    for (i in 1:p)
      rev[i] <- (batx$units[,i] %*% mdiff)^2/(1+ev[i])+log(ev[i]+1/ev[i]+2)
    units <- batx$units[,(order(-rev))]
    proj <- batx$proj[,(order(-rev))]
  }    
  list(ev=ev[order(-rev)], rev=rev[order(-rev)], 
       units=units, proj=proj)
}

# discriminant plot
# bw: black/white
plotcluster <- function(x, clvecd, clnum=NULL,
                        method=ifelse(is.null(clnum),"dc","awc"),bw=FALSE,
                        ignorepoints=FALSE, ignorenum=0, pointsbyclvecd=TRUE,
                        xlab=NULL,
                        ylab=NULL, pch=NULL, col=NULL, ...){
  symvec <- as.integer(clvecd)
  if (min(symvec)==0)
    csymbols <- c(sapply(c(0:9),toString),
      intToUtf8(c(97:122,65:90),multiple=TRUE))
  else
    csymbols <- c(sapply(c(1:9,0),toString),
      intToUtf8(c(97:122,65:90),multiple=TRUE))    
  symvec <- symvec-min(symvec)+1
  if (pointsbyclvecd){
    colsel <- c("black","red","green","blue")
    if (max(symvec>4))
        colsel <- c(colsel,colors()[1+sample(656,max(symvec)-4)])
  }
  if (dim(as.matrix(x))[2]==1){
    clvec <- as.integer(clvecd)
    if (is.null(pch)){
      if (pointsbyclvecd){
        pch <- csymbols[symvec]
        pch[clvecd==ignorenum] <- "N"
      }
      else{
        pch <- symvec
        pch[clvecd==ignorenum] <- 78
      }
    }
    if (is.null(col)){
      col <- if (bw) 1
             else{
               if (pointsbyclvecd)
                 col <- colsel[symvec]
               else
                 col <- symvec
             }
    }
    plot(x,as.integer(clvecd),pch=pch,col=col)    
  }
  else{
    asym <- any(method==c("bc","vbc","adc","awc","arc","anc"))
    if (asym){
      if (is.null(clnum))
        stop("Asymmetric methods require specification of clnum")
      clvec <- as.integer(as.integer(clvecd)==as.integer(clnum))
    }
    else
      clvec <- as.integer(clvecd)
#    print(clvecd)
    cx <- discrproj(x, clvecd, method, clnum, ignorepoints=ignorepoints,
                    ignorenum=ignorenum, ...)$proj
    if (is.null(xlab))
      xlab <- paste(method,"1")
    if (is.null(ylab))
      ylab <- paste(method,"2")
    if (is.null(pch)){
      if (pointsbyclvecd){
        pch <- csymbols[symvec]
        if (ignorepoints)
          pch[clvecd==ignorenum] <- "N"
      }
      else{
        if (asym){
          pch <- 1+clvec
          if (ignorepoints)
            pch[clvecd==ignorenum] <- 78
        }
        else{
          pch <- symvec
          if (ignorepoints)
            pch[clvecd==ignorenum] <- 78
        }
      }
    }
    if (is.null(col)){
      col <- if (bw) 1
             else{
               if (pointsbyclvecd)
                 col <- colsel[symvec]
               else{
                 if (asym)
                   col <- 1+clvec
                 else
                   col <- symvec
               }
             }
    }
    plot(cx, xlab=xlab, ylab=ylab, pch=pch, col=col, ...)
  }
}

discrproj <- function(x, clvecd, method="dc", clnum=NULL,
                      ignorepoints=FALSE, ignorenum=0, ...){
  x <- as.matrix(x)
  if (ignorepoints){
    includei <- clvecd!=ignorenum
    xx <- x[includei,]
    xi <- x[!includei,]
    clvec <- clvecd[includei]
  }
  else{
    xx <- x
    clvec <- clvecd
  }
  result <- switch(method,
                   dc=discrcoord(xx, clvec, ...),
                   bc=batcoord(xx, clvec, clnum),
                   vbc=batcoord(xx, clvec, clnum, dom="var"),
                   mvdc=mvdcoord(xx, clvec, clnum, ...),
                   adc=adcoord(xx, clvec, clnum),
                   awc=awcoord(xx, clvec, clnum, ...),
                   arc=awcoord(xx, clvec, clnum, method="mcd", ...),
                   nc=ncoord(xx, clvec, ...),
                   wnc=ncoord(xx, clvec, weighted=TRUE, ...),
                   anc=ancoord(xx, clvec, clnum, ...))
  if (ignorepoints){
    xip <- xi %*% result$units
    xproj <- x
    xproj[includei,] <- result$proj
    xproj[!includei,] <- xip
    result$proj <- xproj
  }
  result
}
                   
                   







