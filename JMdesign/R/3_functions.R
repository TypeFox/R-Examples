powerLongSurv <- function(N,
                          nevents,
                          tmedian,
                          meantf,
                          p,
                          t,
                          SigmaTheta,
                          sigmae_2,
                          ordtraj,
                          beta = 0, 
                          alpha = 0.05,
                          tol = 1.5e-8){
 

  errorTesting(N=N,
               nevents=nevents,
               tmedian=tmedian,
               meantf=meantf,
               p=p,
               t=t,
               SigmaTheta=SigmaTheta,
               sigmae_2=sigmae_2,
               ordtraj=ordtraj,
               beta=beta, 
               alpha=alpha,
               tol=tol)

 
  # Auxiliary computing functions.
  # *****************************

  ##Calculate integral(x^i eta exp(-eta x) dx) from 0 to a (0<=a<Infinity), 
  ##i=0,1,...,k.
  ##Returns a (k+1) vector.
  integrals <- function(k,a,eta){
                 if( k >= 1 ) {
                   return((1.0 - cumsum(c(1,cumprod((eta*a)/1:k)))*exp(-eta*a))*
                           c(1,cumprod(1:k/eta)) )
                 } else if( k == 0 ) {
                   return( c(1.0 - exp(-eta*a)) )
                 }
               }

  ## Given a matrix A, calculate the sums on diagonals parallel to 
  ## secondary diagonal. 
  diag2Sums <- function(A){
                 A <- as.matrix(A)
                 m <- dim(A)
                 dsums <- numeric(m[1]+m[2]-1)
                 for(j in 1:m[2]) { 
                   dsums[j:(m[1]+j-1)] <- dsums[j:(m[1]+j-1)] + A[,j]
                 }
                 return(dsums)
               }

  f <- function(BSigma, Rn, sigmae_2){
         n <- nrow(Rn)
         tRn <- t(Rn)
         vn <- diag(n)*sigmae_2 + Rn %*% BSigma %*% tRn
         wn <- solve(vn)
         S <- BSigma %*% tRn %*% wn %*% Rn %*% BSigma
         return(S)
       }


  # Main Calculations 
  # *****************
  title <- "Joint Modeling of Longitudinal and Survival Data"
  tnames <- c("Constant", "Linear", "Quadratic", "Cubic", 
              paste(ordtraj, "-th Order", sep="") )
  subtitle <- paste("Power Calculation for Unknown Sigma - ",
                   tnames[min(ordtraj+1,5)], " Trajectory", sep="")

  o <- order(t)
  t <- t[o]
  p <- p[o]
  tt <- c(0,t)

  k <- ordtraj
  BSigma <- SigmaTheta[1:(k+1),1:(k+1)]
    
  # The weighted average of Sigma_theta_hat
  n <- length(tt)
  R <- matrix(1, nrow=n, ncol=1)
  if( k >= 1 ){
    R0 <- matrix(tt, nrow=n, ncol=k, byrow=FALSE)
    R <- cbind(R, R0^rep(1:k,each=n))
  }

  X <- matrix(0, nrow=k+1, ncol=k+1)
  for(i in 2:n){
    Ri <- as.matrix(R[1:i,])
    X <- X + p[i-1] * f(BSigma, Ri, sigmae_2)
  }

  eta <- log(2.0)/tmedian
  censr <- 100.0 * (1.0-nevents/N)
  et <- integrals(k=2*k, a=meantf, eta=eta)
  et[1] <- 1.0
  et[-1] <- (N/nevents)*et[-1]
  zalpha <- qnorm(p={1 - alpha})
  if(beta < 0.0) zalpha <- -zalpha
  mus <- sqrt(diag2Sums(A=X) %*% et) * beta
  prepower <- sqrt(nevents)*mus - qnorm(p={1-alpha/2})
  power <- pnorm(q=prepower)
  if(beta < 0.0)  power <- 1.0 - power

  # output
  # ******
  pls <- new("powerLongSurv",
             title=title,
             subtitle=subtitle,
             t=t,
             p=p,
             N=as.integer(N),
             nevents=as.integer(nevents),
             censr=censr,
             tmedian=tmedian, 
             meantf=meantf,
             SigmaTheta=as.matrix(SigmaTheta), 
             ordtraj=as.integer(ordtraj),
             BSigma=as.matrix(BSigma),
             beta=beta,
             alpha=alpha,
             power=as.numeric(power))

  return(pls)

} ## End of "powerLongSurv()"

