IsoGenem <- function(x, y){

  y <- as.matrix(y)
  if(dim(y)[[1]]==1){
    IsoGene1(x = x, y = y)
  } else {
    
    ordx <- order(x)
    x <- x[ordx]
    y <- y[, ordx]
  
    unx <- unique(x)
    
    ydf <- as.data.frame(t(y))
    y.m <- do.call("cbind", unclass(by(ydf, x, colMeans)))
    
    y.m.tot <- matrix(rep(rowMeans(y), length(x)), ncol = length(x))
    
#    y.is.u <- t(apply(y.m, 1, function(x) isoreg(unx, x)$yf))
#    y.is.d <- t(apply(y.m, 1, function(x) rev(isoreg(rev(unx), x)$yf)))
    n.p <- table(x)
    n.g <- length(n.p) 
    
   y.is.u <- t(apply(y.m, 1, function(x) pava(x, w = n.p)) )
#   y.is.d <- t(apply(y.m, 1, function(x) rev(pava(rev(x), w = rev(n.p))))) 
   y.is.d <- t(apply(y.m, 1, function(x) pava(x, w = n.p, decreasing=TRUE))) 
     
    rep.iso.d <- y.is.d[, rep(1:length(n.p),n.p)]
    rep.iso.u <- y.is.u[, rep(1:length(n.p),n.p)]
    
    y.m.all <- y.m[, rep(1:length(n.p), n.p)]
  
    SST0 <- rowSums((y - rowMeans(y))^2)
  
    SSIS.u1 <- rowSums((rep.iso.u - y)^2)
    SSIS.d1 <- rowSums((rep.iso.d - y)^2)
  
    SST <- rowSums((y - y.m.all)^2)
  
    direction=NULL
  
    direction <- ifelse(SSIS.u1 <= SSIS.d1, "u", "d")
             
    lambda1.up <- SSIS.u1 / SST0
    Esquare.up <- 1 - lambda1.up
    iso.u <- y.is.u
           
    #w.up <- (y.is.u[,n.g] - y.m[,1]) / sqrt(2*SST/(sum(n.p)-n.g)/(n.g-1))
    #w.c.up <- (y.is.u[,n.g] - y.is.u[,1]) / sqrt(2*SST/(sum(n.p)-n.g)/(n.g-1))

    w.up <- (y.is.u[,n.g]-y.m[,1])/sqrt(rowSums((y-y.m.all)^2)/(sum(n.p)-n.g)*(1/n.p[1]+1/n.p[n.g]))   
    w.c.up <- (y.is.u[,n.g]-y.is.u[,1])/sqrt(rowSums((y-y.m.all)^2)/(sum(n.p)-n.g)*(1/n.p[1]+1/n.p[n.g]))   
  

    m.up <- (y.is.u[,n.g] - y.is.u[,1]) / sqrt(SSIS.u1/(sum(n.p)-n.g))
    i.up <- (y.is.u[,n.g] - y.is.u[,1]) /
      sqrt(SSIS.u1/(sum(n.p) - apply(y.is.u, 1, function(x) length(unique(x)))))
  
    lambda1.dn <- SSIS.d1 / SST0
    Esquare.dn <- 1 - lambda1.dn
    iso.u <- y.is.d
  
    n.pSum <- sum(n.p)
    #w.dn <- (y.is.d[,n.g]-y.m[,1])/ sqrt(2*SST / (n.pSum-n.g)/(n.g-1))
    #w.c.dn <- (y.is.d[,n.g]-y.is.d[,1]) / sqrt(2*SST/(n.pSum-n.g)/(n.g-1))
   
     w.dn <- (y.is.d[,n.g]-y.m[,1])/sqrt(rowSums((y-y.m.all)^2)/(sum(n.p)-n.g)*(1/n.p[1]+1/n.p[n.g]))  
     w.c.dn <- (y.is.d[,n.g]-y.is.d[,1])/sqrt(rowSums((y-y.m.all)^2)/(sum(n.p)-n.g)*(1/n.p[1]+1/n.p[n.g]))   
  
   
    m.dn <- (y.is.d[,n.g]-y.is.d[,1]) / sqrt(SSIS.d1/(n.pSum-n.g))
    i.dn <- (y.is.d[,n.g]-y.is.d[,1]) /
      sqrt(SSIS.d1/(sum(n.p) - apply(y.is.d, 1, function(x) length(unique(x)))))
    
    res <-  list(E2.up = Esquare.up,
                 Williams.up = as.numeric(w.up),
                 Marcus.up = as.numeric(w.c.up),
                 M.up = as.numeric(m.up),
                 ModM.up = as.numeric(i.up),
                 E2.dn = Esquare.dn,
                 Williams.dn = as.numeric(w.dn),
                 Marcus.dn = as.numeric(w.c.dn),
                 M.dn = as.numeric(m.dn),
                 ModM.dn = as.numeric(i.dn),
                 direction = direction)
    return(res)
  }
}
