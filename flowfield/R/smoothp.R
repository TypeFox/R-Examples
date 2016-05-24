smoothp <- function (t,y,x,d) {
  #*********************************************************************************************
  #  This function determines the smoothing parameter lambda for the penalized spline 
  #  regression.  The process has 2 steps.  First use Wand's asymptotic approximate to obtain an 
  #  initial guess and then use this an input to the optimize function which minimizes the cross 
  #  validation procedure.
  #
  #    Refererences:  1. D. Ruppert, M. P. Wand and R. J. Carroll, Semiparametric Regression. 
  #                      New York, NY: Cambidge University Press, 2003.
  #
  #                   2. M. P. Wand, On the optimal amount of smoothing in penalized spline 
  #                      regression, Biometrika 86 (1999), 936–940. 
  #
  #                   3. Hutchinson, Michael F., and F. R. De Hoog. "Smoothing noisy data with 
  #                      spline functions." Numerische Mathematik 47.1 (1985): 99-106. 
  #
  #    Input: t - time series observation times 
  #           y - time series response values
  #           x - the design matrix for the penalized spline regression (PSR)
  #           d - the diagonal matrix used for smoothing the PSR
  #           
  #    Output: The PSR smoothing parameter lambda
  #
  #*********************************************************************************************
  u <- length(t)
  knots <- floor(u/10)
  space <- (t[u] - t[1])/knots 
  
  pvals <- rep(0,u)  # Vector to hold the first smoothed approximation using quadratic 
  # polynomials. 
  
  svals <- rep(0,u)  # Vector to hold the variance of the first smoothed approximation.
  vpos <- 1          # Pointer to hold vector position  
  
  #**********************************************************************************************
  # Chunk the time series observations into knot intervals and fit quadratic polynomials
  # the data.
  #**********************************************************************************************
  for (i in 1:(knots+1)){
    j <- 1
    tchunk <- t[vpos+j-1]
    ypoly <- y[vpos+j-1]
    if (i==1){limit <- t[1] + space/2}
    if (i>1) {limit <- limit + space} 
    while((vpos+j <= length(t)) && (t[vpos+j] <= limit)) {
      tchunk <- c(tchunk,t[vpos+j])
      ypoly <- c(ypoly,y[vpos+j])
      j <- j+1
    }
    vpos <- vpos + length(tchunk)
    n <- length(tchunk)
    degree <- 2
    if(n < 3) {
      write("Warning: Data is too sparse in some sections, forecast results unreliable",file="")
    }
    model <- lm(ypoly ~ poly(tchunk,degree))
    tchunk <- data.frame(tchunk)
    if (i==1){
      pvals <- predict(model,tchunk)
      svals <- sum(resid(model)^2)/(n-(degree+1))
    }
    if(i>1) {
      pvals <- c(pvals,predict(model,tchunk))
      svals <- c(svals,sum(resid(model)^2)/(n-(degree+1)))
    }
  }
  
  pvals <- data.frame(pvals)
  svals <- data.frame(svals)
  
  sigma2 <- mean(svals[,])
  
  num <- sigma2*sum(diag(solve(t(x)%*%x+10^(-10)*d,d)))
  denom <- norm(x%*%solve(t(x)%*%x+10^(-10)*d,d)%*%solve(t(x)%*%x+10^(-10)*d,t(x))%*%pvals[,1],type="2")^2
  
  guess <- sqrt(num/denom) # Initial guess for lambda based on Wand's asymptotic approximation.
  library(stats)
  # Determine lambda by minimizing the generalized cross validation criteria (Hutchinson and de Hoog)  
  lambda <- optimize(GCV,interval=0:(3*guess),lower=0,upper=3*guess,tol=0.01,y=y,x=x,d=d) 
  
  lambda <- lambda$minimum
  
  # Returns the guess if optimize didn't converge
  if ((lambda == 0) || (lambda == 3*guess)) {lambda <- guess}
  
  return(lambda)
}

