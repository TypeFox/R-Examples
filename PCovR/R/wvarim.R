wvarim <-
function(F1,nrs=20){
  n <- nrow(F1)
  m <- ncol(F1)
  conv <- 1e-06
  
  #normalize rows of F1
  Q <- eigen(F1 %*% t(F1))$vectors
  D <- eigen(F1 %*% t(F1))$values
  F <- Q[,1:m] %*% sqrt(diag(D[1:m]))
  #reflect rows with negative first factor loadings
  tau <- diag(sign(F[,1]))
  H <- diag(diag(F %*% t(F))^(-1/2))
  A <- tau %*% H %*% F
  #symmetry about 1st axis of A
  c <- sqrt(1/m)
  
  acos_a <- acos(A[,1])*360/(2*pi)
  acos_c <- matrix(acos(c)*360/(2*pi),n,1)
  W1 <- A[,1]<c
  M1 <- cos(((acos_a-acos_c)/(matrix(90,n,1)-acos_c)*90)*2*pi/360)^2+matrix(0.001,n,1)
  W2 <- A[,1]>=c
  M2 <- cos(((acos_a-acos_c)/(acos_c)*90)*2*pi/360)^2+matrix(0.001,n,1)
  W <- W1*M1+W2*M2 
  A <- Diagonal(n,W) %*% A
  
  r <- ncol(A)
  for (nos in 1:(nrs+1)){
    if (nos==1){
      T <- diag(r)
    } else {
      T <- orth(matrix(rnorm(r*r),r,r))
    }
    B <- A %*% T
    f <- SUM(scale(B*B, center = TRUE, scale = FALSE))$ssq
    fold <- f-2*conv*f
    if (f==0){
      fold <- -conv
    }
    iter <- 0
    while(f-fold > f*conv) {
      fold <- f
      iter <- iter+1
      for (i in 1:(r-1)){
        for (j in (i+1):r){
          x <- B[,i]
          y <- B[,j]
          xx <- T[,i]
          yy <- T[,j]
          u <- x^2-y^2
          v <- 2*x*y
          u <- scale(u, center = TRUE, scale = FALSE)
          v <- scale(v, center = TRUE, scale = FALSE)
          a <- 2*sum(u*v)
          b <- sum(u^2)-sum(v^2)
          c <- sqrt(a^2+b^2)
          if (a >= 0){
            Sign <- 1
          } else {
            Sign=-1
          }
          if (c<10e-12){
            stop("rotation does not make sense for single factor models.")
          } else {
            vvv <- -Sign*sqrt(((b+c)/(2*c)))
            sinus <- sqrt(.5-.5*vvv)
            cosinus <- sqrt(.5+.5*vvv)
          }
          v <- cosinus*x-sinus*y
          w <- cosinus*y+sinus*x
          vv <- cosinus*xx-sinus*yy
          ww <- cosinus*yy+sinus*xx
          if (vvv >= 0){
            B[,i] <- v
            B[,j] <- w
            T[,i] <- vv
            T[,j] <- ww
          } else {
            B[,i] <- w
            B[,j] <- v
            T[,i] <- ww
            T[,j] <- vv            
          }
        }
      }
      f <- SUM(scale(B*B, center = TRUE, scale = FALSE))$ssq
    }
    if (nos==1 || f>fr){
      fr <- f
      ir <- iter
      Tr <- T
    }
  }
  I <- diag(r)
  A <- F %*% Tr
  for (i in 1:r){
    if (sum(A[,i])<0){
      I[i,i]=-1
    }
  }
  
  Tr <- Tr %*% I
  A <- F %*% Tr
  Tr <- solve(t(F1) %*% F1) %*% t(F1) %*% A
  loadings <- F1 %*% Tr
  rot <- list(Th=Tr,loadings=loadings,W=W,fr=fr,ir=ir)
  return(rot)
}
