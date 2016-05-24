delta.error.loop2r <- function (g) {
if (g$method=="harmonic2") {
  z <- g$fit[[1]]
  rss <- sum(z$residuals^2)
  p <- z$rank
  p1 <- 1L:p
  resvar1 <- rss/z$df.residual  
  R1 <- chol2inv(z$qr$qr[p1, p1, drop = FALSE])
  se1 <- sqrt(diag(R1) * resvar1)
  
  
  z2 <- g$fit[[2]]
  rss <- sum(z2$residuals^2)
  p <- z2$rank
  p1 <- 1L:p
  resvar2 <- rss/z2$df.residual  
  R2 <- chol2inv(z2$qr$qr[p1, p1, drop = FALSE])
  se2 <- sqrt(diag(R2) * resvar2)
  
  cov.Ta<- R1*resvar1
  cov.Tb<- R2*resvar2
  
  cov.matrix <- cbind(rbind(cov.Ta,matrix(0,4,3)),rbind(matrix(0,3,4),cov.Tb))
  retention.aboveSE <- cov.matrix[5,5]+cov.matrix[7,7]+2*cov.matrix[5,7]
  b.xSE<-deltamethod(~sqrt(x2^2+x3^2),c(z$coefficients,z2$coefficients), cov.matrix)
  phase.angleSE<-deltamethod(~atan(x3/x2),c(z$coefficients,z2$coefficients), cov.matrix)

  if (g$values["n"]==1) splitSE <-   deltamethod(~atan(x6/sqrt(x2^2+x3^2)),c(z$coefficients,z2$coefficients), cov.matrix)*180/pi
else splitSE <- NA
  hysteresis.y.belowSE <- deltamethod(~x5/x6,c(z$coefficients,z2$coefficients), cov.matrix)
    hysteresis.y.aboveSE <- deltamethod(~(x5+x7)/x6,c(z$coefficients,z2$coefficients), cov.matrix)

  m <- g$values["m"]
  form <- sprintf("~1/sqrt(1+(x6/x5)^(2/%f))",m)
  hysteresis.x.belowSE<-deltamethod(as.formula(form),c(z$coefficients,z2$coefficients), cov.matrix)
  form <- sprintf("~1/sqrt(1+(x6/(x5+x7))^(2/%f))",m)
  hysteresis.x.aboveSE<-deltamethod(as.formula(form),c(z$coefficients,z2$coefficients), cov.matrix)
  
  coercionform <- sprintf("~sqrt(x2^2+x3^2)/sqrt(1+(x6/x5)^(2/%f))",m)
  coercion.belowSE<-deltamethod(as.formula(coercionform),c(z$coefficients,z2$coefficients), cov.matrix)
    coercionformabove <- sprintf("~sqrt(x2^2+x3^2)/sqrt(1+(x6/(x5+x7))^(2/%f))",m)
  coercion.aboveSE<-deltamethod(as.formula(coercionform),c(z$coefficients,z2$coefficients), cov.matrix)
  lag.belowSE <- deltamethod(~atan(x5/x6),c(z$coefficients,z2$coefficients), cov.matrix)*g$period/(2*pi)
    lag.aboveSE <- deltamethod(~atan((x5+x7)/x6),c(z$coefficients,z2$coefficients), cov.matrix)*g$period/(2*pi)

  leftpart <- (0.5/(beta((m + 3)/2, (m + 3)/2) * (m + 2)) + 1/beta((m + 
          1)/2, (m + 1)/2) - 1/beta((m + 3)/2, (m - 1)/2))/(2^m) * pi
  form <- sprintf("~(%f * (x5+x7/2)* atan(x3/x2))/2",leftpart)
  areaSE <-   deltamethod(as.formula(form),c(z$coefficients,z2$coefficients), cov.matrix)

  SEs<- list("n"=NA,"m"=NA,"b.x"=b.xSE,"b.y"=se2[3],"phase.angle"=phase.angleSE,
             "cx"=se1[1],"cy"=se2[1],"retention.above"=retention.aboveSE,"retention.below"=se2[2],"coercion.above"=coercion.aboveSE,"coercion.below"=coercion.belowSE,
             "area"=areaSE,"lag.above"=lag.aboveSE,"lag.below"=lag.belowSE, "split.angle"=splitSE,"hysteresis.x.above"=hysteresis.x.aboveSE,
             "hysteresis.x.below"=hysteresis.x.belowSE,"hysteresis.y.above"=hysteresis.y.aboveSE,"hysteresis.y.below"=hysteresis.y.belowSE)
  }
  else {
    n <- length(g$x) 
    coefs <- g$fit$par[(n+1):(n+8)]
    vmat <- (as.vector(2*crossprod(g$residuals)/(n-8))*solve(g$fit$hessian))
    vmat2 <- vmat[(n+1):(n+8),(n+1):(n+8)]
    SEm <- deltamethod(~exp(x5),coefs,vmat2)
    SEn <- deltamethod(~exp(x6),coefs,vmat2)
    coercion.aboveSE <- deltamethod(~x3/sqrt(1+(x4/x7)^(2/x5)),coefs,vmat2)
    coercion.belowSE <- deltamethod(~x3/sqrt(1+(x4/x8)^(2/x5)),coefs,vmat2)
    hysteresis.x.aboveSE <- deltamethod(~1/sqrt(1+(x4/x7)^(2/x5)),coefs,vmat2)
    hysteresis.x.belowSE <- deltamethod(~1/sqrt(1+(x4/x8)^(2/x5)),coefs,vmat2)
    areaSE <- NA
    splitSE <- NA
    hysteresis.y.aboveSE <- deltamethod(~x7/x4,coefs,vmat2)
    hysteresis.y.belowSE <- deltamethod(~x8/x4,coefs,vmat2)
    lag.aboveSE <- deltamethod(~atan(x7/x4),coefs,vmat2)*g$period/(pi*2)
    lag.belowSE <- deltamethod(~atan(x8/x4),coefs,vmat2)*g$period/(pi*2)
    SEs<- list("n"=SEn,"m"=SEm,"b.x"=sqrt(vmat2[3,3]),"b.y"=sqrt(vmat2[4,4]),"phase.angle"=sqrt(vmat[1,1]),
               "cx"=sqrt(vmat2[1,1]),"cy"=sqrt(vmat2[2,2]),"retention.above"=sqrt(vmat2[7,7]),"retention.below"=sqrt(vmat2[8,8]),
               "coercion.above"=coercion.aboveSE,"coercion.below"=coercion.belowSE,"area"=areaSE,"lag.above"=lag.aboveSE,
               "lag.below"=lag.belowSE,"split.angle"=splitSE,"hysteresis.x.above"=hysteresis.x.aboveSE,
               "hysteresis.x.below"=hysteresis.x.belowSE,"hysteresis.y.above"=hysteresis.y.aboveSE,"hysteresis.y.below"=hysteresis.y.belowSE)
    
  }
  SEs
}
