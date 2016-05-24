CADFpvalues <- function(t0, rho2=0.5, type=c("trend", "drift", "none"))
{
  # This procedure computes the p-values of the CADF test developed by Hansen (1995)
  # t0:        sample statistic
  # rho2:      value of the parameter rho^2
  # type:     CADF model type. No constant "nc", Constant "c", Constant plus trend "ct"
  # The procedure is described in Costantini, Lupi & Popp (2007)
  # Citation:
  # @TECHREPORT{,
  # author = {Costantini, Mauro and Lupi, Claudio and Popp, Stephan},
  # title = {A Panel-{CADF} Test for Unit Roots},
  # institution = {University of Molise},
  # year = {2007},
  # type = {Economics \& Statistics Discussion Paper},
  # number = {39/07},
  # timestamp = {2008.10.15},
  # url = {http://econpapers.repec.org/paper/molecsdps/esdp07039.htm}}

  type <- match.arg(type)

  switch(type,
         "trend" = coeffs <- coeffs_ct,
         "drift" = coeffs <- coeffs_c,
         "none"  = coeffs <- coeffs_nc)

  # the first column of coefs are the probabilities, the other columns are beta_0, ..., beta_3
  # of Costantini, Lupi & Popp eqn (13).
    
  L <- dim(coeffs)[1]

  # compute the fitted quantiles for the given value of rho^2
  fitted.q <- coeffs[,2] + coeffs[,3]*rho2 + coeffs[,4]*rho2^2 + coeffs[,5]*rho2^3

  # find the position of the fitted quantile that is closest to the sample statistic
  difference <- abs(fitted.q - t0)
  position   <- which(difference==min(difference))
  if (length(position)>1) position <- position[1]
  
  # interpolate locally using eq (10) in Costantini, Lupi & Popp (2007) using l observations (l must be an integer odd number)
  l <- 11
  if ( (position > ((l-1)/2)) & (position < (L - (l-1)/2)) ) range <- (position - (l-1)/2):(position + (l-1)/2)
  if (position <= ((l-1)/2))                                 range <- 1:l
  if (position >= (L-(l-1)/2))                               range <- (L-(l-1)/2+1):L
  prob <- coeffs[range,1] 
  prob <- qnorm(prob)
  x1 <- fitted.q[range]
  x2 <- x1^2
  x3 <- x1^3
  local.interp <- lm(prob~x1+x2+x3)
  model.summary <- summary(local.interp)
  gamma  <- model.summary$coefficients
  vt0 <- c(1, t0, t0^2, t0^3)
  vt0 <- vt0[1:dim(gamma)[1]]
  p.value <- vt0%*%gamma[,1]
  p.value <- as.vector(pnorm(p.value))
  return(p.value)
}
