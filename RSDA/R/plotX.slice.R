plotX.slice <-
function(xx1,yy1,xx2,yy2,vv,vvars,kk) {
  radio <- sqrt(xx1^2+yy1^2)
  a <- xx1
  b <- yy1
  cx1 <- 0
  cy1 <- 0
  if(radio <= 1.0) {
    cx1 <- xx1
    cy1 <- yy1
    if(xx1>0) {
      text(xx1+0.15, yy1,vvars[kk],col=vv[kk])
    }
    else {
      text(xx1-0.15, yy1,vvars[kk],col=vv[kk])
    }
  }
  else {
    na <- -(a/sqrt(a^2 + b^2))
    nb <- -(b/sqrt(a^2 + b^2))
    sig1 <- a*na
    sig2 <- b*nb
    if((sig1 < 0) & (sig2 < 0)) {
      na <- -na
      nb <- -nb
    }
    cx1 <- na
    cy1 <- nb
    if(na>0) {
      text(na+0.15, nb,vvars[kk],col=vv[kk])
    }
    else {
      text(na-0.15, nb,vvars[kk],col=vv[kk])
    }
  }
  radio <- sqrt(xx2^2+yy2^2)	
  a <- xx2
  b <- yy2
  cx2 <- 0
  cy2 <- 0
  if(radio <= 1.0) {
    cx2 <- xx2
    cy2 <- yy2
  }
  else {
    na <- -(a/sqrt(a^2 + b^2))
    nb <- -(b/sqrt(a^2 + b^2))
    sig1 <- a*na
    sig2 <- b*nb
    if((sig1 < 0) & (sig2 < 0)) {
      na <- -na
      nb <- -nb
    }
    cx2 <- na
    cy2 <- nb
  }
  fxx <- c(0,cx1,cx2,0)
  fyy <- c(0,cy1,cy2,0)
  polygon(fxx, fyy, col=vv[kk], border = vv[kk])
}
