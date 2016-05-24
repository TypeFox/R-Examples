direct <- function(x,y,w=NULL) {
  n <- length(x)
  if (!is.null(w))
    W <- diag(c(sqrt(w)))
      else W <- diag(n)
  D1 <- matrix(ncol=3,nrow=n)
  D1[,1] <- x^2
  D1[,2] <- x*y
  D1[,3] <- y^2
  D1 <- W%*%D1
  D2 <- W%*%cbind(x,y,rep(1,n))
  
  S1 <- crossprod(D1)
  S2 <- crossprod(D1,D2)
  S3 <- crossprod(D2)
  
  c0<-numeric(9)
  c0[3] <- -2
  c0[5] <- 1
  c0[7]<- -2
  C1 <- matrix( c0, ncol=3)
  
  Tmatrix <- - solve(S3,t(S2))
  M <- S1 + S2 %*% Tmatrix
  M <- solve(C1,M)
  
  gev <- eigen(M)
  gevec <- gev$vectors
  cond <- 4*gevec[1,]*gevec[3,]-gevec[2,]^2
  a <- numeric(6)
  a[1:3] <- gevec[,cond > 0]
  a[4:6] <- Tmatrix %*% a[1:3]
  theta <- atan2(a[2],a[1]-a[3])/2
  while(theta<0){theta<-pi/2+theta} 
  
  cx <- -(2*a[3]*a[4]-a[2]*a[5])/(4*a[1]*a[3]-a[2]*a[2])
  cy <- -(2*a[1]*a[5]-a[2]*a[4])/(4*a[1]*a[3]-a[2]*a[2])
  
  major <- 1/sqrt((a[1]*cos(theta)*cos(theta) + a[2]*cos(theta)*sin(theta) + a[3]*sin(theta)*sin(theta)) / (a[1]*cx*cx + a[2]*cx*cy + a[3]*cy*cy - a[6]))
  minor <- 1/sqrt((a[1]*sin(theta)*sin(theta) - a[2]*cos(theta)*sin(theta) + a[3]*cos(theta)*cos(theta)) / (a[1]*cx*cx + a[2]*cx*cy + a[3]*cy*cy - a[6]))
  
  
  semi.major <- abs(major) 
  semi.minor <- abs(minor)
  
  if (semi.minor > semi.major){
    semi.minor <- semi.major; semi.major <- abs(minor); theta <- theta +pi/2;
  }
  rotated.angle <- theta*180/pi
  res <- cbind(D1,D2)%*%a
  sigma2 <- sum(res^2)/(n-5)
  sigmaa <- sigma2*ginv(crossprod(cbind(D1,D2)))
    names(a) <- c("x2","xy","y2","x","y","int")
  list("residuals"=res,"fit"=list("a"=a,"eigen"=gev,"sigma2"=sigma2,"vcov"=sigmaa),vals=c("cx"=cx,
                "cy"=cy,"theta"=theta,"semi.major"=semi.major,"semi.minor"=semi.minor,"rotated.angle"=rotated.angle))
}
