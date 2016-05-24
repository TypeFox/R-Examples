##################################################################
#
#
# FUNCTIONS TO CALCULATE UNCERTAINTY IN BACKGOUND ESTIMATION
#
#


##################################################################
# get.hessian(data, knots.x, knots.y, Gr, r, p.bkg)
#
# computes Hessian matrix of the target function
# for lazy calculation get.hess.numerically is recommended 
# as more stable function
get.hess <- function(data, knots.x, knots.y, Gr=NA, r=seq(0, 2, 0.01), p.bkg=0.5){
  x <- data$x
  y <- data$y-data$SB
  sigma <- data$sigma
  lambda <- data$lambda
  Phi <- basisMatrix(x=x, knots.x=knots.x)
  bkg <- Phi %*% t(t(knots.y))
  
  # 1. Prior 
#  cat("Calculating prior hessian... \n")
  D <- DMatrix(knots.x=knots.x)$matrix  
  cDc <- as.vector(t(knots.y) %*% D %*% t(t(knots.y)))
  E <- length(knots.y)

  bkg.pp <- D %*% t(t(knots.y)) # second derivative of bkg
  hess.prior <- E/2 * (2*(D / cDc) - (4 * bkg.pp %*% t(bkg.pp)) / (cDc ^ 2) )

  # 2. Likelihood:
 # cat("Calculating likelihood hessian...")  	
  deviation <- y - bkg
  deviation.norm <- deviation/sigma
  grad.f <- as.vector(-deviation/sigma^2)
  hess.f <- as.vector(-1/sigma^2)
  
  funct.f <- (log(p.bkg) - 0.5 * log(2 * pi) - log(sigma) - 0.5 * deviation.norm ^ 2)
  f <- list(funct=funct.f, grad=grad.f, hess=hess.f)

   #######
  rho <- sigma/lambda
  z <- (y - bkg) / lambda
  qq <- z / rho - rho
  funct.h <- log(1 - p.bkg) - log(lambda) + pnorm(log.p=TRUE, q=qq) - z  + 0.5 * (rho ^ 2)
 
  gamma.q <- exp(-0.5 * qq ^ 2 - pnorm(q=qq, log.p=TRUE)) / (rho * sqrt(2 * pi))
  grad.h <- as.vector(-(1 - gamma.q) / lambda) 
  
  hess.h <- as.vector(-gamma.q * (gamma.q + qq / rho) / (lambda ^ 2))
  h <- list(funct=funct.h, grad=grad.h, hess=hess.h)
  
   #######    
  f.fract <- as.vector(1 / (1 + exp(h$funct - f$funct)))
  h.fract <- 1 - f.fract
  grad.g <- hess.g <- NA

  grad.g <- (f.fract * f$grad + (1 - f.fract) * h$grad) 
  f.contr <- f.fract*f$hess + f.fract*f$grad*f$grad
  h.contr <- h.fract*h$hess + h.fract*h$grad*h$grad
  
  hess.gg <- -(f.contr + h.contr - (grad.g)^2)
  
  Phi.prime <- hess.gg*Phi
  hess.g <- t(Phi.prime) %*% Phi
 
  # 3. Gr=-4*Pi*rho*r restriction
  hess.gr <- 0
  if(!is.na(Gr[1])){
  #  cat("Calculating Gr hessian...")   
	if(is.na(Gr$type1))
      hess.gr.r1 <- 0
    else if(Gr$type1=="gaussianNoise")
      hess.gr.r1 <- logLikelihoodGrGauss(y=data$y-data$SB, knots.y=knots.y, alpha=1, Phi=Phi, bkg.r=Gr$bkg.r, 
	                                 sigma.r=Gr$sigma.r, matrix.FT=Gr$matrix.FT1, Hessian=TRUE)$hess
	else if(Gr$type1=="correlatedNoise")
	  hess.gr.r1 <- logLikelihoodGrCorr(knots.y=knots.y, Phi=Phi, bkg.r=Gr$bkg.r, 
	                                KG.inv=Gr$KG.inv, matrix.FT=Gr$matrix.FT1, Hessian=TRUE)$hess
									
	if(is.na(Gr$type2))
      hess.gr.r2 <- 0
	else if(Gr$type2=="secondDeriv")
	  hess.gr.r2 <- logPriorBkgRSmooth(bkg.r=Gr$matrix.FT2 %*% bkg, D=Gr$D, Hessian=TRUE, Phi=Phi, 
	                                   matrix.FT=Gr$matrix.FT2, knots.y=knots.y)$hess
	else if(Gr$type2=="gaussianProcess")
      hess.gr.r2 <- logPriorBkgRGP(bkg.r=Gr$matrix.FT2 %*% bkg, covMatrix=Gr$covMatrix, Hessian=TRUE)$hess								
	hess.gr <- hess.gr.r1 + hess.gr.r2														
#	cat(" done!\n")								
  }
  
  # 4. Summing up
  hess <- hess.prior + hess.gr + hess.g
  #hess <- hess.g

  hess.inv <- solve(hess)
   
  # 5. Converting Hessian into Q-space from a spline-space  
  Phi <- basisMatrix(x=data$x, knots.x=knots.x)
  H <- Phi%*%hess.inv%*%t(Phi)   # hessian in Q-space (inverted...)
                                 # that is, covariance matrix
  H <- H + diag((data$sigma)^2, length(data$x))
  
  cov.diag <- diag(H)
  cov.diag[which(cov.diag<0)]<-0
  stdev <- sqrt(cov.diag)

  # 6. Converting Hessian into r-space
  MFT <- sineFT.matrix(Q=data$x, r=r)
  cov.r <- MFT %*% H %*% t(MFT)
  cov.diag.r <- diag(cov.r)
  cov.diag.r[which(cov.diag.r<0)]<-0
  stdev.r <- sqrt(cov.diag.r)

#  return(list(stdev=stdev, stdev.r=stdev.r, hess=hess, cov.matrix=hess.inv, cov.matrix.r=cov.r, hess.gg=hess.gg, H=H))
  return(list(stdev=stdev, stdev.r=stdev.r, hess=hess, cov.matrix=hess.inv, cov.matrix.r=cov.r, hess.gg=hess.gg))


}

##################################################################
#get.hess.numerically(data, knots.x, knots.y, Gr, r, p.bkg, h1)
#
# computes a list with elements
#   stdev:          estimated standard deviations for a reconstructed signal.
#   stdev.r:        estimated standard deviations for a reconstructed signal in r-space.
#   hess:           Hessian matrix for a target function.
#   cov.matrix:     covariance matrix, i.e. the inverse of the Hessian.
#   cov.matrix.r:   covariance matrix in r-space.
get.hess.numerically <- function(data, knots.x, knots.y, Gr=NA, r=seq(0, 2, 0.01), p.bkg=0.5, h=1e-4){
  n <- length(knots.y)
  incr <- function(cc, h, i=0, j=0){
	if(i!=0)
	  cc[i] <- cc[i]+h
	if(j!=0)
	  cc[j] <- cc[j]+h
	return(cc)
  }
  
  # 1. Computing Hessian matrix
  hess <- matrix(0, nrow=n, ncol=n)
  for(i in 1:n){
  cat("knot.i = ", i, " of ", n,"\n")
    for(j in 1:n){ 
	  a1 <- logPosterior(data=data, alpha=1, knots.x=knots.x, knots.y=incr(cc=knots.y,h=h,i=i,j=j), Gr=Gr, p.bkg=p.bkg) 
	  a2 <- logPosterior(data=data, alpha=1, knots.x=knots.x, knots.y=incr(cc=knots.y,h=h,i=i, j=0), Gr=Gr, p.bkg=p.bkg) 
	  a3 <- logPosterior(data=data, alpha=1, knots.x=knots.x, knots.y=incr(cc=knots.y,h=h,i=0, j=j), Gr=Gr, p.bkg=p.bkg) 
	  a4 <- logPosterior(data=data, alpha=1, knots.x=knots.x, knots.y=knots.y, Gr=Gr, p.bkg=p.bkg) 
	  hess[i,j] <- (a1-a2-a3+a4)/h^2
	}  
  } 
  # 2. Computing inverse 
  # U <- regularized.cholesky(hess)
  # hess.inv <- chol2inv(U)
  hess.inv <- solve(hess)
   
  # 3. Converting Hessian into Q-space from a spline-space  
  Phi <- basisMatrix(x=data$x, knots.x=knots.x)
  H <- Phi%*%hess.inv%*%t(Phi)   # hessian in Q-space (inverted...)
                                 # that is, covariance matrix
  H <- H + diag((data$sigma)^2, length(data$x))
  
  cov.diag <- diag(H)
  cov.diag[which(cov.diag<0)]<-0
  stdev <- sqrt(cov.diag)

  # 4. Converting Hessian into r-space
  MFT <- sineFT.matrix(Q=data$x, r=r)
  cov.r <- MFT %*% H %*% t(MFT)
  cov.diag.r <- diag(cov.r)
  cov.diag.r[which(cov.diag.r<0)]<-0
  stdev.r <- sqrt(cov.diag.r)

  return(list(stdev=stdev, stdev.r=stdev.r, hess=hess, cov.matrix=hess.inv, cov.matrix.r=cov.r))
}


##############################################################################################
# grad.descent(data, knots.x, knots.y, Gr, p.bkg, N)
#
# Gradient descent method to find a local minimum for a psi(c)
grad.descent <- function(data, knots.x, knots.y, Gr=NA, p.bkg=0.5, eps=1e-3, N=10000){
  x.p <- knots.y
  x.pp <- x.p
  lambda <- c(abs(0.0001/get.deriv(data=data, knots.x, x.p, Gr=Gr, p.bkg=0.5)))
  N <- round(N, digits=-2)
  if(N==0) N <- 100

  cat("iterations will stop once convergence is reached \n")
  cat("indicating % of itermax... \n")
  for(j in 1:(N/100)){
    cat("...",(j-1)/N*100*100, "% done \n")
    for(i in 1:100){
	  grad.f <- as.vector(get.deriv(data=data, knots.x=knots.x, knots.y=x.p, Gr=Gr, p.bkg=p.bkg))
#	grad.f <- as.vector(get.deriv.numerically(data=data, knots.x=knots.x, knots.y=x.p, Gr=Gr, p.bkg=p.bkg))
      x.f <- x.p - lambda*grad.f  
      x.p <- x.f
    }
	if(max(abs(x.pp-x.f))<eps){
	  cat("\n convergence reached! \n")
	  break
	}
	x.pp <- x.f
  }
  if(j == (N/100)){ 
    cat("convergence not reached! \n")
	  return(knots.y)
  }
  else
    return(x.f)
}



#####################################################################################
# get.deriv.numerically(data, knots.x, knots.y, Gr, p.bkg,  h)
#
# Function to calculate derivate of psi(c) over c. 
# More stable than get.deriv
get.deriv.numerically <- function(data, knots.x, knots.y, Gr=NA, p.bkg=0.5,  h=1e-6){
  n <- length(knots.y)
  incr <- function(cc, h, i=0){
    cc[i] <- cc[i]+h
	cc	
  }
  der <- 0
  for(i in 1:n){
    a1 <- logPosterior(data=data, alpha=1, knots.x=knots.x, knots.y=incr(cc=knots.y,h=h,i=i), Gr=Gr, p.bkg=p.bkg) 
    a2 <- logPosterior(data=data, alpha=1, knots.x=knots.x, knots.y=knots.y, Gr=Gr, p.bkg=p.bkg) 
    der[i] <- (a1-a2)/h
  } 
 return(der)
}


#####################################################################################
# get.deriv(data, knots.x, knots.y, Gr, p.bkg,)
#
# Function to calculate derivate of psi(c) over c. 
# Less stable than get.deriv.numerically
get.deriv <- function(data, knots.x, knots.y, Gr=NA, p.bkg=0.5){
  x <- data$x
  y <- data$y-data$SB
  sigma <- data$sigma
  lambda <- data$lambda
  
  Phi <- basisMatrix(x=x, knots.x=knots.x)
  bkg <- Phi %*% t(t(knots.y))
  
 # 1. Prior 
 # cat("Calculating prior hessian... \n")
  D <- DMatrix(knots.x=knots.x)$matrix  
  cDc <- as.vector(t(knots.y) %*% D %*% t(t(knots.y)))
  E <- length(knots.y)


  # 2. Likelihood:
  deviation <- y - bkg
  norm.dev <- deviation/sigma
  grad.f <- as.vector(deviation/sigma^2)*Phi

  funct.f <- (log(p.bkg) - 0.5 * log(2 * pi) - log(sigma) - 0.5 * norm.dev ^ 2)
  f <- list(funct=funct.f, grad=grad.f)
  
  
  #######
  rho <- sigma/lambda
  z <- (y - bkg) / lambda
  qq <- z / rho - rho
  funct.h <- (log(1 - p.bkg) - log(lambda) + pnorm(log.p=TRUE, q=qq) - z  + 0.5 * (rho ^ 2))

  gamma.q <- exp(-0.5 * qq ^ 2 - pnorm(q=qq, log.p=TRUE)) / (rho * sqrt(2 * pi))
  grad.h <- as.vector((1 - gamma.q) / lambda) * Phi
  h <- list(funct=funct.h, grad=grad.h)
  
  #######    
  f.fract <- as.vector(1 / (1 + exp(h$funct - f$funct)))
  grad.g <- hess.g <- NA
  grad.g <- f.fract * f$grad + (1 - f.fract) * h$grad

  # 4. Summing up
  grad.g <- f.fract * f$grad + (1 - f.fract) * h$grad
  grad <- -colSums(grad.g) + (t(knots.y) %*% D) * E / cDc
	
	
  ########
  if(!is.na(Gr[1])){
    matrix.FT <- Gr$matrix.FT1
	sigma.r <- Gr$sigma.r
  
    Mprime <- matrix.FT%*%Phi / (sqrt(2)*sigma.r)
	b.prime <- Gr$bkg.r / (sqrt(2)*sigma.r)
    grad.gr <- 2*t(Mprime) %*% Mprime %*%knots.y- 2*t(Mprime)%*%b.prime
    grad <- as.vector(grad) + as.vector(grad.gr)
  }
	
  return(grad)
}
 


# author: Eric Cai
# http://www.r-bloggers.com/scripts-and-functions-using-r-to-implement-the-golden
#                           -section-search-method-for-numerical-optimization/
golden.search = function(data, lower.bound=-.01, upper.bound=.01, tolerance=1e-6, knots.x, knots.y, Gr, p.bkg, grad.f){

  f <- function(lambda){
    logPosterior(data=data, alpha=1, knots.x=knots.x, knots.y=(knots.y - lambda*grad.f), Gr=Gr, p.bkg=p.bkg)
  }
 
  golden.ratio = 2/(sqrt(5) + 1)
  x1 = upper.bound - golden.ratio*(upper.bound - lower.bound)
  x2 = lower.bound + golden.ratio*(upper.bound - lower.bound)
  f1 = f(x1)
  f2 = f(x2)
  iteration = 0
  while(abs(upper.bound - lower.bound) > tolerance){
    iteration = iteration + 1
    if (f2 > f1){
      upper.bound = x2
      x2 = x1
      f2 = f1
      x1 = upper.bound - golden.ratio*(upper.bound - lower.bound)
      f1 = f(x1)
    } 
    else{
      lower.bound = x1
      x1 = x2
      f1 = f2
      x2 = lower.bound + golden.ratio*(upper.bound - lower.bound)
      f2 = f(x2)
    }
  }
  estimated.minimizer = (lower.bound + upper.bound)/2
  estimated.minimizer  
}

#####################################################################################
# regularized.cholesky(Matr, eps.max, eps.min, numTries) 
#
# Regularized Cholesky matrix decomposition
regularized.cholesky <- function(Matr, eps.max=1e-2, eps.min=1e-20, numTries=17) {
  baseVal <- min(diag(Matr))
  U <- try(chol(Matr), silent=T)
  epsilon <- eps.min
  I <- diag(nrow(Matr))
  while (!is.null(attr(U, "class")) && attr(U, "class") == "try-error" && epsilon <= eps.max) {
    U <- try(chol(Matr + baseVal * epsilon * I), silent=T)
    epsilon <- epsilon * 10
  }
  if (epsilon >= eps.max) stop("We just couldn't Cholesky-decompose this matrix\n")
  return (U)
}


row.outer.product <- function(Phi) {
  # Computes the (3d) array consisting of the outer product of each row of Phi
  # with itself.
  #
  # Args:
  #   Phi:  A (R x C) numeric matrix (usually a spline matrix, where the
  #      columns correspond to the knots, and the rows correspond to evaluation
  #      points).
  #
  # Returns:
  #   A 3-index numeric array, whose dimensions are c(C, C, R), such that the
  #   matrix at [*, *, i] is the outer product of Phi[i, ] with itself.
  rows <- nrow(Phi)
  cols <- ncol(Phi)
  ppt <- array(apply(Phi, 1, function(x) x %*% t(x)), dim=c(cols, cols, rows))
  return (ppt)
}

#########################################################
# covMatrixSE(x, sig, l, noiseFactor)
#
# returns:   covariance matrix for a squared-exponential cov function
# arguments 
#   x:       datapoints
#   sig:     vertical scale for variations
#   l:       horizontal scale for variations: correlation length
#            either a sigle value or vector of length 'length(x)'
covMatrixSE <- function(x, sig=0.05, l=0.1){
  N <- length(x)
  covX <- matrix(nrow=N, ncol=N)
  if(length(l)==1){
    covX <- outer(X=x, Y=x, FUN=function(x, y) {
     sig^2*exp( -0.5*((x - y) / l)^2 )
    })
  }
  else{ # make it faster: l[i]^2+l[j]^2 calculate only one time
        # also use symmetric properties
    for(i in 1:N)
      for(j in 1:N)
	    covX[i,j] <- sig^2*sqrt( 2*l[i]*l[j]/(l[i]^2+l[j]^2) )*   
		          exp( -(x[i]-x[j])^2/(l[i]^2+l[j]^2) )  
  }
 # covX <- covX - sig^2 + diag(sig^2, N)
  factor <- 1
  while (det(covX)==0) {
    covX <- covX*1e1
    factor <- factor*1e1
  }
  
  list(cov=covX, factor=factor)
}

#########################################################
# covMatrix.DI(covMatrix)
#
# returns:      determinant and inverse of the covariance matrix
# arguments 
#   covMatrix:  covariance matrix
covMatrix.DI <- function(covMatrix){

  U <- regularized.cholesky(covMatrix)  # change to carefulChol
  covMatrix.inv <- chol2inv(U) 
  covMatrix.det <- det(U)
  
  list(inv=covMatrix.inv, det=covMatrix.det)
}