#'@title Likelihood Estimation for liv
#'@description Computes the log-likelihood function. Only two groups are considered, since as presented in Ebbes et al (2005) this gives good, unbiased results.
#'@param theta - a vector of initial values for the parameters of the model to be supplied to the optimization algorithm.
#'@param y - a vector or matrix containing the dependent variable.
#'@param P -  a vector with the endogeneous variable or a matrix of dimention n X 2, where each column contains an endogeneous variable
#'@return returns the value of the negative log-likelihood.
#'@keywords endogenous
#'@keywords latent
#'@keywords instruments
#'@keywords internal
#'@author  adapted by Raluca Gui from the code provided by Professor Ebbes during a workshop at Univ. of Zurich in April 2015.
#'@references Ebbes, P., Wedel,M., Boeckenholt, U., and Steerneman, A. G. M. (2005). 'Solving and testing for regressor-error
#'(in)dependence when no instruments
#'
logL <- function(theta,y,P){
  # parameters
  b00 <- theta[1]
  a1 <- theta[2]

  # means of the groups
  pi1 <- theta[3]
  pi2 <- theta[4]

    # Reparametricise to avoid neg. variances or abs(correlation) >1
    # var e = a^2
    #	var v = b^2 + c^2
    #	cov ev = ab
    #	rho ev = ab / (a + sqrt(b^2 + c^2))


    Sigma_temp <- matrix(0,nrow=2, ncol=2)
    Sigma_temp[1,1] <- theta[5]
    Sigma_temp[1,2] <- 0
    Sigma_temp[2,1] <- theta[6]
    Sigma_temp[2,2] <- theta[7]

    Sigma  <-  Sigma_temp %*% t(Sigma_temp)

    s2e <- Sigma[1,1]
    s2v <- Sigma[2,2]
    sev <- Sigma[1,2]

    pt <- exp(theta[8])   # Probability for group 1

    pt <- pt / (1+pt)

  # Build varcov matrix reduced form
  varcov  <-  matrix(0,2,2)
  varcov[1,1] <- a1*a1*s2v+2*a1*sev+s2e
  varcov[2,1] <- a1*s2v+sev
  varcov[1,2] <- varcov[2,1]
  varcov[2,2] <- s2v


  # Group 1 contribution
  mu1 <- matrix(nrow=2,ncol=1)
  mu1[1,1] <- b00+a1*pi1
  mu1[2,1] <- pi1

  pdf1 <- mvtnorm::dmvnorm(cbind(y,P), mean=mu1, sigma=varcov)


  # Group 2 contribution
  mu2 <- matrix(nrow=2,ncol=1)
  mu2[1,1] <- b00+a1*pi2
  mu2[2,1] <- pi2

  pdf2 <- mvtnorm::dmvnorm(cbind(y,P), mean=mu2, sigma=varcov)


  logll <-  sum(log(pt*pdf1 + (1-pt)*pdf2))

  return(-1*logll)


}#End likelihood function


