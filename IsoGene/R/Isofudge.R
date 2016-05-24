#' calculate the fudge factor
#' 
#' @param x.dif mean difference 
#' @param si2  sd  
xfudge <- function(x.dif, si2) {
  xs.alpha <- quantile(si2, (0:20) / 20)
  xs.perc <- quantile(si2, (0:100) / 100)
  xs.perc.gr <- as.numeric(cut(si2, breaks=xs.perc, include.lowest=TRUE))
  xs.mad <- matrix(0, length(xs.alpha), length(xs.perc)-1)
  for (ii in 1:100){
    k0 <- xs.perc.gr == ii
    xs.mad[,ii] <- apply(x.dif[k0]/outer(si2[k0],xs.alpha,"+"),2,mad, constant = 1)
  }
  xcv <- function(x) sqrt(var(x)) 
  xs.cv.sd <- apply(xs.mad, 1, xcv)
  xs.cv.m <- rowMeans(xs.mad)
  xs.cv <- xs.cv.sd/xs.cv.m
  xfudge <- xs.alpha[sort.list(xs.cv)[1]]
  return(xfudge)
}


Isofudge <- function(x,y){
  
  y <- as.matrix(y)
  
  if(dim(y)[[1]]==1){
    warning("There is only one gene in the data set.")
  } else {  
    ordx <- order(x) # compute once
    x <- x[ordx]
    y <- y[,ordx]  # reverse order (sorted two times)
    
    unx <- unique(x) # compute once
    
    ydf <- as.data.frame(t(y))
    y.m <- do.call("cbind", unclass(by(ydf, x, colMeans)))
    
    y.m.tot <- matrix(rep(rowMeans(y), length(x)), ncol = length(x))

    n.p <- table(x)
    n.g <- length(n.p) 
    
  y.is.u <- t(apply(y.m, 1, function(x) pava(x, w = n.p)) )
  y.is.d <- t(apply(y.m, 1, function(x) rev(pava(rev(x), w = rev(n.p))))) 
 
    
   # y.is.u <- t(apply(y.m, 1, function(x) isoreg(unx, x)$yf))
  #  y.is.d <- t(apply(y.m, 1, function(x) rev(isoreg(rev(unx), x)$yf)))
    

    
    ###################################################
    
    rep.iso.d <- y.is.d[, rep(1:length(n.p),n.p)]
    rep.iso.u <- y.is.u[, rep(1:length(n.p),n.p)]
    
    y.m.all <- y.m[, rep(1:length(n.p), n.p)]
    
    ########################################################
    
    SST0 <- rowSums((y - rowMeans(y))^2)
    
    SSIS.u1 <- rowSums((rep.iso.u - y)^2)
    SSIS.d1 <- rowSums((rep.iso.d - y)^2)
    
    SST <- rowSums((y - y.m.all)^2)
    
    direction <- NULL
    
    direction <- ifelse(SSIS.u1 <= SSIS.d1, "u", "d")
    
    iso.u <- y.m
    SSIS.dir <- NULL
    
    iso.u[direction=="u",] <- y.is.u[direction=="u",]
    iso.u[direction=="d",] <- y.is.d[direction=="d",]
    
    SSIS.dir[direction=="u"] <- SSIS.u1[direction=="u"]
    SSIS.dir[direction=="d"] <- SSIS.d1[direction=="d"]
    
    E2.dif <- sqrt(SST0-SSIS.dir)
    E2.si2 <- sqrt(SST0)
    
    W.dif <- iso.u[,n.g] - y.m[,1]
    W.si2 <- sqrt(SST/(sum(n.p)-n.g)*(1/n.p[1] + 1/n.p[n.g])) 
    
    W.C.dif <- iso.u[,n.g] - iso.u[,1]
    W.C.si2 <- sqrt(SST/(sum(n.p)-n.g)*(1/n.p[1] + 1/n.p[n.g])) 
    
    M.dif <- iso.u[,n.g] - iso.u[,1]
    M.si2 <- sqrt(SSIS.dir/(sum(n.p)-n.g)) 
    
    I.dif <- iso.u[,n.g] - iso.u[,1]
    I.si2 <- sqrt(SSIS.dir/(sum(n.p) - apply(iso.u, 1, function(x) length(unique(x))))) 
    
    fudge.E2 <- xfudge(E2.dif,E2.si2)
    fudge.Williams <- xfudge(W.dif, W.si2)
    fudge.Marcus <- xfudge(W.C.dif,W.C.si2)
    fudge.M <- xfudge(M.dif,M.si2)
    fudge.I <- xfudge(I.dif, I.si2)  
    
    return(as.numeric(c(fudge.E2, fudge.Williams, fudge.Marcus, fudge.M, fudge.I)))
  }
}
