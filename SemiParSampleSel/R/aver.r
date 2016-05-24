aver <- function(x, sig.lev = 0.05, sw = NULL, univariate = FALSE, delta = TRUE, n.sim = 100){

bd <- x$bd
sigma <- x$sigma
nu <- x$nu
X2sg <- x$X2s
X2.d2 <- x$X2.d2
X3.d2 <- x$X3.d2
X4.d2 <- x$X4.d2
X3 <- x$X3
X4 <- x$X4

if (is.null(X3)==TRUE) X3 <- matrix(1, dim(X2sg)[1], 1); X3.d2 <- 1
if (is.null(X4)==TRUE) X4 <- matrix(1, dim(X2sg)[1], 1); X4.d2 <- 1


if(univariate==TRUE) { 
     if(x$margins[2] %in% c("N","G") ) eta2 <- X2sg%*%coef(x$gam2) else eta2 <- X2sg%*%x$gamlss.fit$fit$argument[1:x$X2.d2]  
}

if(univariate==FALSE) eta2 <- x$eta2



if(is.null(sw)) sw <- rep(1,length(eta2)) 
if(x$margins[2] == "N") core <- apply( c(eta2)*X2sg,      2, weighted.mean,  w = sw) 
if(x$margins[2] %in% c("G", "P", "GEOM", "YULE", "NB", "PIG", "NBII", "WARING", "ZIP2", "D", "S") ) 
  { core <- apply( c(exp(eta2))*X2sg, 2, weighted.mean,  w = sw) }
if(x$margins[2] %in% c("BI", "BB") )
  { core <- apply( c( plogis(eta2)*bd )*X2sg, 2, weighted.mean,  w = sw) }
if(x$margins[2] %in% c("LG") )
  { core <- apply( c( ((-(log(1-plogis(eta2)))^(-1))*plogis(eta2))/(1-plogis(eta2)) )*X2sg, 2, weighted.mean,  w = sw) }
if(x$margins[2] %in% c("ZIBI") )
  { core <- apply( c( (1-sigma)*plogis(eta2)*bd )*X2sg, 2, weighted.mean,  w = sw) }
if(x$margins[2] %in% c("ZABI") )
  { core <- apply( c( ((1-sigma)*plogis(eta2)*bd)/(1-(1-plogis(eta2))^bd) )*X2sg, 2, weighted.mean,  w = sw) }
if(x$margins[2] %in% c("ZIBB") )
  { core <- apply( c( (1-nu)*plogis(eta2)*bd )*X2sg, 2, weighted.mean,  w = sw) }
if(x$margins[2] %in% c("ZABB") )
  { core <- apply( c( ((1-nu)*plogis(eta2)*bd)/(1-dBB(0, mu=plogis(eta2), sigma=sigma, bd=bd)) )*X2sg, 2, weighted.mean,  w = sw) }
if(x$margins[2] %in% c("ZINBI", "ZIPIG") )
  { core <- apply( c( (1-nu)*exp(eta2) )*X2sg, 2, weighted.mean,  w = sw) }
if(x$margins[2] %in% c("ZANBI") )
  { core <- apply( c( ((1-nu)*exp(eta2))/(1-(1+sigma*exp(eta2))^(-1/sigma)) )*X2sg, 2, weighted.mean,  w = sw) }
if(x$margins[2] %in% c("ZIP") )
  { core <- apply( c( (1-sigma)*exp(eta2) )*X2sg, 2, weighted.mean,  w = sw) }
if(x$margins[2] %in% c("ZAP") )
  { core <- apply( c( ((1-sigma)*exp(eta2))/(1-exp(-exp(eta2))) )*X2sg, 2, weighted.mean,  w = sw) }
if(x$margins[2] %in% c("ZALG") )
  { core <- apply( c( ((1-sigma)*(-(log(1-plogis(eta2)))^(-1))*plogis(eta2))/(1-plogis(eta2)) )*X2sg, 2, weighted.mean,  w = sw) }



if(univariate==FALSE){

if(!is.null(x$X3) ){
  if(x$margins[2] %in% c("N","G","NB", "PIG", "BB", "NBII", "WARING", "ZABI", "ZIBI", "ZALG", "ZAP", "ZIP", "ZIP2") )  zs <- rep(0, x$X3.d2 + x$X4.d2)  
  if(x$margins[2] %in% c("P", "BI", "GEOM", "LG", "YULE") )                                                            zs <- rep(0, x$X3.d2) 
  if(x$margins[2] %in% c("D", "S", "ZIBB", "ZABB", "ZANBI", "ZINBI", "ZIPIG") )                                        zs <- rep(0, x$X3.d2 + x$X4.d2 + x$X5.d2)  
                   }  

if(is.null(x$X3) ){
  if(x$margins[2] %in% c("N","G","NB", "PIG", "BB", "NBII", "WARING", "ZABI", "ZIBI", "ZALG", "ZAP", "ZIP", "ZIP2") )  zs <- c(0, 0)  
  if(x$margins[2] %in% c("P", "BI", "GEOM", "LG", "YULE") )                                                            zs <- c(0) 
  if(x$margins[2] %in% c("D", "S", "ZIBB", "ZABB", "ZANBI", "ZINBI", "ZIPIG") )                                        zs <- c(0, 0, 0)  
                  }          
}


if(univariate==TRUE){

if(!is.null(x$X3) ){
  if(x$margins[2] %in% c("N","G","NB", "PIG", "BB", "NBII", "WARING", "ZABI", "ZIBI", "ZALG", "ZAP", "ZIP", "ZIP2") )  zs <- rep(0, x$X3.d2)  
  if(x$margins[2] %in% c("P", "BI", "GEOM", "LG", "YULE") )                                                            zs <- NULL 
  if(x$margins[2] %in% c("D", "S", "ZIBB", "ZABB", "ZANBI", "ZINBI", "ZIPIG") )                                        zs <- rep(0, x$X3.d2 + x$X4.d2)  
                   }  

if(is.null(x$X3) ){
  if(x$margins[2] %in% c("N","G","NB", "PIG", "BB", "NBII", "WARING", "ZABI", "ZIBI", "ZALG", "ZAP", "ZIP", "ZIP2") )  zs <- c(0)  
  if(x$margins[2] %in% c("P", "BI", "GEOM", "LG", "YULE") )                                                            zs <- NULL
  if(x$margins[2] %in% c("D", "S", "ZIBB", "ZABB", "ZANBI", "ZINBI", "ZIPIG") )                                        zs <- c(0, 0)  
                  }          
}


if(univariate==FALSE)                                    G <- c( rep(0, x$X1.d2), core , zs)  
if(univariate==TRUE &&   x$margins[2] %in% c("N","G") )  G <- c( core) 
if(univariate==TRUE && !(x$margins[2] %in% c("N","G")))  G <- c( core, zs ) 
  
  
  
if(x$margins[2] == "N") wm <- weighted.mean(eta2, w=sw)
if(x$margins[2] %in% c("G", "P", "GEOM", "YULE", "NB", "PIG", "NBII", "WARING", "ZIP2", "D", "S") )
  { wm <- weighted.mean(exp(eta2), w=sw) }
if(x$margins[2] %in% c("BI", "BB") )
  { wm <- weighted.mean(plogis(eta2)*bd, w=sw) }
if(x$margins[2] %in% c("LG") )
  { wm <- weighted.mean(((-(log(1-plogis(eta2)))^(-1))*plogis(eta2))/(1-plogis(eta2)), w=sw) }
if(x$margins[2] %in% c("ZIBI") )
  { wm <- weighted.mean( (1-sigma)*plogis(eta2)*bd, w=sw) }
if(x$margins[2] %in% c("ZABI") )
  { wm <- weighted.mean( ((1-sigma)*plogis(eta2)*bd)/(1-(1-plogis(eta2))^bd), w=sw) }
if(x$margins[2] %in% c("ZIBB") )
  { wm <- weighted.mean( (1-nu)*plogis(eta2)*bd, w=sw) }  
if(x$margins[2] %in% c("ZABB") )
  { wm <- weighted.mean(((1-nu)*plogis(eta2)*bd)/(1-dBB(0, mu=plogis(eta2), sigma=sigma, bd=bd)), w=sw) }
if(x$margins[2] %in% c("ZINBI", "ZIPIG") )
  { wm <- weighted.mean( (1-nu)*exp(eta2), w=sw) }
if(x$margins[2] %in% c("ZANBI") )
  { wm <- weighted.mean( ((1-nu)*exp(eta2))/(1-(1+sigma*exp(eta2))^(-1/sigma)), w=sw) }
if(x$margins[2] %in% c("ZIP") )
  { wm <- weighted.mean( (1-sigma)*exp(eta2), w=sw) }
if(x$margins[2] %in% c("ZAP") )
  { wm <- weighted.mean( ((1-sigma)*exp(eta2))/(1-exp(-exp(eta2))), w=sw) }
if(x$margins[2] %in% c("ZALG") )
  { wm <- weighted.mean( ((1-sigma)*(-(log(1-plogis(eta2)))^(-1))*plogis(eta2))/(1-plogis(eta2)), w=sw) }




if(univariate==FALSE) Vv <- x$Vb

if(univariate==TRUE){
	if(x$margins[2] %in% c("N","G") ) Vv <- x$gam2$Vp  else {
	   hess <- x$gamlss$fit$hessian                                                                                                     
	   He.eig <- eigen(hess, symmetric=TRUE)
	   if(min(He.eig$values) < sqrt(.Machine$double.eps)) He.eig$values[which(He.eig$values < sqrt(.Machine$double.eps))] <- 0.0000001
	   if(length(He.eig$values) == 1) {
        Vv <- He.eig$vectors%*%tcrossprod(1/He.eig$values,He.eig$vectors)
      } else {
        Vv <- He.eig$vectors%*%tcrossprod(diag(1/He.eig$values),He.eig$vectors)
      }
	}
}


rMVN <- function (n, mean, sigma) {
    L <- mroot(sigma)
    lL <- ncol(L)
    t(mean + L %*% matrix(rnorm(lL * n), lL, n))
}
    
if (delta == TRUE) {
    sv <- sqrt( t(G)%*%Vv%*%G ) 
    qz <- qnorm(sig.lev/2, lower.tail = FALSE)
    lb <- wm - qz*sv 
    ub <- wm + qz*sv 
    }
if (delta == FALSE) {
    if (univariate == FALSE) 
      coefm <- x$coefficients
    if(univariate==TRUE) { 
      coefm <- x$gamlss.fit$fit$argument
      if(x$margins[2] %in% c("N") ) {coefm <- coef(x$gam2)} 
      if (x$margins[2] %in% c("G") ) {coefm <- coef(x$gam2)}
    }
    bs <- rMVN(n.sim, mean = coefm, sigma = Vv)
    bs.mu <- bs[ , 1:x$X2.d2]
    if (univariate == FALSE)
      bs.mu <- bs[, x$X1.d2 + (1:x$X2.d2)]
if(x$margins[2] == "N") p2s <- X2sg %*% t(bs.mu)
if(x$margins[2] %in% c("G", "P", "GEOM", "YULE", "NB", "PIG", "NBII", "WARING", "ZIP2", "D", "S") )
  { p2s <- exp(X2sg %*% t(bs.mu)) }
if(x$margins[2] %in% c("BI", "BB") )
  { p2s <- plogis(X2sg %*% t(bs.mu))*bd }
if(x$margins[2] %in% c("LG") )
  { p2s <- ((-(log(1-plogis(X2sg %*% t(bs.mu))))^(-1))*plogis(X2sg %*% t(bs.mu)))/(1-plogis(X2sg %*% t(bs.mu))) }
if(x$margins[2] %in% c("ZIBI") ){ 
  if (univariate == FALSE) {bs.sigma <- bs[, x$X1.d2 + X2.d2 + (1:X3.d2)]} else { bs.sigma <- bs[, X2.d2 + (1:X3.d2)] }
  post.sigma <- plogis(X3 %*% t(bs.sigma))
  p2s <- (1-post.sigma)*plogis(X2sg %*% t(bs.mu))*bd 
  }
if(x$margins[2] %in% c("ZABI") ){ 
  if (univariate == FALSE) {bs.sigma <- bs[, x$X1.d2 + X2.d2 + (1:X3.d2)]} else { bs.sigma <- bs[, X2.d2 + (1:X3.d2)] }
  post.sigma <- plogis(X3 %*% t(bs.sigma))
  p2s <-  ((1-post.sigma)*plogis(X2sg %*% t(bs.mu))*bd)/(1-(1-plogis(X2sg %*% t(bs.mu)))^bd) 
}
if(x$margins[2] %in% c("ZIBB") ){ 
  if (univariate == FALSE) {bs.nu <- bs[, x$X1.d2 + X2.d2 + X3.d2 + (1:X4.d2)]} else { bs.nu <- bs[, X2.d2 + X3.d2 + (1:X4.d2)] }
  post.nu <- plogis(X4 %*% t(bs.nu))
  p2s <-  (1-post.nu)*plogis(X2sg %*% t(bs.mu))*bd 
  }  
if(x$margins[2] %in% c("ZABB") ) { 
  if (univariate == FALSE) {bs.sigma <- bs[, x$X1.d2 + X2.d2 + (1:X3.d2)]} else { bs.sigma <- bs[, X2.d2 + (1:X3.d2)] }
  post.sigma <- plogis(X3 %*% t(bs.sigma))
  if (univariate == FALSE) {bs.nu <- bs[, x$X1.d2 + X2.d2 + X3.d2 + (1:X4.d2)]} else { bs.nu <- bs[, X2.d2 + X3.d2 + (1:X4.d2)] }
  post.nu <- plogis(X4 %*% t(bs.nu))
  p2s <- ((1-post.nu)*plogis(X2sg %*% t(bs.mu))*bd)/(1-dBB(0, mu=plogis(X2sg %*% t(bs.mu)), sigma=post.sigma, bd=bd)) 
  }
if(x$margins[2] %in% c("ZINBI", "ZIPIG") ){ 
  #sw <- sw[!x$y2==0]  does excluding zeros make a difference?
  #X2sg <- X2sg[!x$y2==0, ]
  #X4 <- X4[!x$y2==0, ]
  if (univariate == FALSE) {bs.nu <- bs[, x$X1.d2 + X2.d2 + X3.d2 + (1:X4.d2)]} else { bs.nu <- bs[, X2.d2 + X3.d2 + (1:X4.d2)] }
  post.nu <- plogis(X4 %*% t(bs.nu))
  p2s <-  (1-post.nu)*exp(X2sg %*% t(bs.mu)) 
  }
if(x$margins[2] %in% c("ZANBI") ){ 
  if (univariate == FALSE) {bs.sigma <- bs[, x$X1.d2 + X2.d2 + (1:X3.d2)]} else { bs.sigma <- bs[, X2.d2 + (1:X3.d2)] }
  post.sigma <- plogis(X3 %*% t(bs.sigma))
  if (univariate == FALSE) {bs.nu <- bs[, x$X1.d2 + X2.d2 + X3.d2 + (1:X4.d2)]} else { bs.nu <- bs[, X2.d2 + X3.d2 + (1:X4.d2)] }
  post.nu <- plogis(X4 %*% t(bs.nu))
  p2s <- ((1-post.nu)*exp(X2sg %*% t(bs.mu)))/(1-(1+post.sigma*exp(X2sg %*% t(bs.mu)))^(-1/post.sigma)) 
  }
if(x$margins[2] %in% c("ZIP") )  { 
  if (univariate == FALSE) {bs.sigma <- bs[, x$X1.d2 + X2.d2 + (1:X3.d2)]} else { bs.sigma <- bs[, X2.d2 + (1:X3.d2)] }
  post.sigma <- plogis(X3 %*% t(bs.sigma))
  p2s <- (1-post.sigma)*exp(X2sg %*% t(bs.mu)) 
  }
if(x$margins[2] %in% c("ZAP") ){ 
  if (univariate == FALSE) {bs.sigma <- bs[, x$X1.d2 + X2.d2 + (1:X3.d2)]} else { bs.sigma <- bs[, X2.d2 + (1:X3.d2)] }
  post.sigma <- plogis(X3 %*% t(bs.sigma))
  p2s <- ((1-post.sigma)*exp(X2sg %*% t(bs.mu)))/(1-exp(-exp(X2sg %*% t(bs.mu)))) 
  }
if(x$margins[2] %in% c("ZALG") ){
  if (univariate == FALSE) {bs.sigma <- bs[, x$X1.d2 + X2.d2 + (1:X3.d2)]} else { bs.sigma <- bs[, X2.d2 + (1:X3.d2)] }
  post.sigma <- plogis(X3 %*% t(bs.sigma))
  p2s <- ((1-post.sigma)*(-(log(1-plogis(X2sg %*% t(bs.mu))))^(-1))*plogis(X2sg %*% t(bs.mu)))/(1-plogis(X2sg %*% t(bs.mu))) 
  }
    wms <- apply(p2s, MARGIN = 2, FUN = weighted.mean, 
      w = sw)
    bb <- quantile(wms, probs = c(sig.lev/2, 1 - sig.lev/2), 
        na.rm = TRUE)
    lb <- bb[1]
    ub <- bb[2]
}


  res <- c(lb, wm, ub)

  out <- list(res = res, sig.lev = sig.lev)
 
  class(out) <- "aver"

  out

}


# importFrom(gamlss.dist,dNBI,pNBI,dDEL,pDEL,dPIG,pPIG,dSI,pSI,tofyDEL2)
