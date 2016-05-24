"logit.hessian" <-
function(x=x,f=f,del=rep(.002,length(x$est)),dapprox=FALSE,nfcn=0) {
  # computes the hessian of f on logit scale
  small <- 1.e-8
  npar <- length(x$est)
  nlm <- 2*npar*npar
  np2 <- 2*npar
  
  xt <- ftrf(x$est, x$low, x$upp)
  f0 <- f(x$est); nfcn <- nfcn + 1
  # cat('compute internal Hessian with function value at:',f0,'\n')

  # *** INITIALIZE NEWTON BY CALCULATION OF A
  #     NLM POINT SIMPLEX GEOMETRY IN N DIM PARAMETER SPACE
  xn <- matrix(rep(xt,nlm),ncol=nlm)
  
  # ***   ON AXES
  # ***	compute delta distance in each direction

  for(i in 1:npar) {
  j.even <- 2*i
  j.odd  <- j.even-1
  xn[i,j.even] <- xn[i,j.even]-del[i]
  xn[i,j.odd]  <- xn[i,j.odd] +del[i]
  }

  if(dapprox==FALSE) {
   # ***   OFF AXIS
   if(npar > 1) {
     mc <- np2+1
     for(i in 2:npar) {
       for(j in 1:(i-1)) {
         xn[i,mc] <- xn[i,mc]+del[i]
         xn[j,mc] <- xn[j,mc]+del[j]
         mc <- mc+1
         xn[i,mc] <- xn[i,mc]-del[i]
         xn[j,mc] <- xn[j,mc]+del[j]
         mc <- mc+1
         xn[i,mc] <- xn[i,mc]-del[i]
         xn[j,mc] <- xn[j,mc]-del[j]
         mc <- mc+1
         xn[i,mc] <- xn[i,mc]+del[i]
         xn[j,mc] <- xn[j,mc]-del[j]
         mc <- mc+1
       }}
   }
   
   f.vec <- numeric(nlm) 
   for(i in 1:nlm) {
     f.vec[i] <- f(btrf(xn[,i], x$low, x$upp))
     nfcn <- nfcn + 1
   }
 } else { # ignore off diagonal elements

   f.vec <- numeric(np2) 
   for(i in 1:np2) {
     f.vec[i] <- f(btrf(xn[,i], x$low, x$upp))
     nfcn <- nfcn + 1
   }
 }
  
  # ***   FIRST AND DIAGONAL SECOND DERIVATIVES
 
  i <- 1:npar; i.even <- 2*i; i.odd <- i.even-1
 
  df <- (f.vec[i.even]-f.vec[i.odd])/2/del
  if(npar > 1) {ddf <- diag((f.vec[i.even]+f.vec[i.odd]-2*f0)/(del**2))/2 }
  else {ddf <- (f.vec[i.even]+f.vec[i.odd]-2*f0)/(del**2)/2 }  
  # print(format(ddf),quote=FALSE)
 
  if(dapprox==FALSE) {
  # ***   SECOND DERIVATIVES
 
  if(npar > 1) {
 
    mc <- np2+1
    for(i in 2:npar) {
      for(j in 1:(i-1)) {
        mc1 <- mc; mc2 <- mc1+1; mc3 <- mc2+1; mc4 <- mc3+1
        ddf[i,j] <- (f.vec[mc1]+f.vec[mc3]-f.vec[mc2]-f.vec[mc4])/del[i]/del[j]/4
        mc <- mc4+1
      }}
    ddf <- ddf+t(ddf)
  }
} else {ddf <- 2 * ddf}
 
  # ***   EIGEN VALUES and MATRIX INVERSION
 
  eig <- eigen(ddf)
  # cat('eigen values: ',format(sort(eig$values)),'\n')
  if(any(eig$values < 0)) {
    warning('hessian not pos. definite')
  }
  if(any(abs(eig$values) < small)) {
    warning('hessian may be singular')
  }
  # *** ADJUSTMENT OF del (feature not implemented) 
  # ddf.diag <- diag(ddf)
  # del[ddf.diag > 0] <- .002/sqrt(ddf.diag[ddf.diag > 0])

  return(list(df=df,ddf=ddf,nfcn=nfcn,eigen=eig$values))
}
