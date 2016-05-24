BoxCox <- function(y, lambda, rescale=TRUE, na.rm=rescale){
##
## 1.  lambda, sign.y, log(abs(y+lambda[2]))
##
  eps <- (1/50) # from MASS::boxcox 
  la <- lambda[1]
  if(length(lambda)<2){
    y0 <- 0
    sgn.y <- sign(y, 1)
    lay <- log(abs(y))
  } else {
    y0 <- lambda[2]  
    sgn.y <- sign(y+lambda[2], 1)
    lay <- log(abs(y+lambda[2]))
  }
##
## 2.  geometric mean
##
  if(is.numeric(rescale)){
    g <- rescale
    rescale <- TRUE 
  } else if(na.rm){
    g <- exp(mean(lay[!is.na(y)]))
  } else g <- exp(mean(lay)) 
##
## 3.  BoxCox
##
  laly <- lambda[1]*lay 
#  bc <- (sgn.y*expm1(laly)+(sgn.y-1))/la
  bc <- (sgn.y*expm1(laly))/la
  bc[sgn.y<0] <- (bc[sgn.y<0]-(2/la))
  smy <- (abs(laly)<eps)
  if(any(smy)){
#   finite Tayor approx to expm1(la*ly)/ly  
# =     
# = sgn.y*ly*(1+(laly/2)*(1+(laly/3)*(1+(laly/4)*(1+...))))
#        + (sgn.y-1)/la}
    lasm <- laly[smy]
    bs <- lasm/8
    for(j in 7:2){
      bs <- (bs+1)*lasm/j
    }
#    bs. <- (sgn.y[smy]*lay[smy]*(bs+1) + (sgn.y[smy]-1)/la) 
    bs. <- sgn.y[smy]*lay[smy]*(bs+1) 
    bs.[sgn.y[smy]<0] <- (bs.[sgn.y[smy]<0] - (2/la))
    bc[smy] <- bs.
  }
##
## 4.  rescale 
##  
  if(rescale){
    bc <- bc/(g^(la-1))
  } 
##
## 5.  attr 
##
  attr(bc, 'lambda') <- lambda
  if(any(sgn.y<1)) attr(bc, 'sign.y') <- sgn.y 
  attr(bc, 'rescale') <- rescale 
  attr(bc, 'GeometricMean') <- g  
  class(bc) <- 'BoxCox'
##
## 6.  Done 
## 
  bc
} 

invBoxCox <- function(z, lambda, sign.y, GeometricMean, rescale){
##
## 1.  lambda
##
  eps <- (1/50) # from MASS::boxcox 
  if(missing(lambda)){
    lambda <- attr(z, 'lambda')
    if(is.null(lambda))
      stop('lambda missing with no default')
  }
  la <- lambda[1]
##
## 2.  sign.y 
##
  if(missing(sign.y)){
    sign.y <- attr(z, 'sign.y')
    if(is.null(sign.y))sign.y <- 1    
  }
##
## 3.  GeometricMean 
##
  if(missing(GeometricMean)){
    GeometricMean <- attr(z, 'GeometricMean')
    if(is.null(GeometricMean))GeometricMean <- 1 
  }
##
## 4.  rescale 
##
  if(missing(rescale)){
    rescale <- attr(z, 'rescale')
    if(is.null(rescale))rescale <- TRUE 
# note:  With GeometricMean=1, 
# the answer is the same with rescale = TRUE or FALSE 
  }
##
## 5.  w = z*GeometricMean^(la-1)
##
  w <- as.vector(z)
  if(rescale){
    w <- w*GeometricMean^(la-1)
  }  
##
## 6.  lay = log(abs(y+y0))
##
  law <- la*w 
  lay <- (log(abs(1+law))/la)
  smlw <- (abs(law)<eps)
  if(any(smlw)){
#  finite Taylor approx to log1p(la*w)/la
#  = w*(1+la*w*(1/2)+la*w*((1/3)+la*w*(1+...)))
    ay <- 1
    lws <- law[smlw]
    for(j in 8:1){
      ay <- lws*((1/j)-ay)
    }
    lay[smlw] <- (w[smlw]*ay)    
  }
##
## 7.  y <- exp(lay)+y0 
##
  ay <- exp(lay)
  if(length(lambda)>1){
    ay <- ay - lambda[2]
  } 
  y <- sign.y*ay 
##
## 8.  Done 
##
  y
}