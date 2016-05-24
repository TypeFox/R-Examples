dirichlet.moments <-
function# Calculate the moments of the Dirichlet distribution given the values of vector alpha
### internal function for sisus
(alpha.in
### internal variable
, priors.precision
### internal variable
, priors.sources
### internal variable
)
{
  ##details<<
  ## interal function for sisus.run()

  # dirEV() by Erik Barry Erhardt 12/26/2006 3:40PM
  # this function is for exploration of priors on the p vector
  # Given a vector of dirichlet parameters alpha, it will return the first couple moments of the p's

  alpha = alpha.in * priors.precision;
  alpha0 = sum(alpha);
  p.E = alpha/alpha0;                                 # Expectations
  p.Var = alpha*(alpha0-alpha)/(alpha0^2*(alpha0+1)); # Variance
  p.SD = sqrt(p.Var);                                 # Standard Deviation

  p.o = paste("           Dirichlet prior on p vector",                                         "\n"); write.out(p.o);
  p.o = paste("           Input priors:    ",           paste(priors.sources,   collapse=", "), "\n"); write.out(p.o);
  p.o = paste("           Input precision: ",           paste(priors.precision, collapse=", "), "\n"); write.out(p.o);
  p.o = paste("           Alphas:   ",                  paste(alpha,            collapse=", "), "\n"); write.out(p.o);
  p.o = paste("           Expected: ",                  paste(p.E,              collapse=", "), "\n"); write.out(p.o);
  p.o = paste("           Variance: ",                  paste(p.Var,            collapse=", "), "\n"); write.out(p.o);
  p.o = paste("           StdDev:   ",                  paste(p.SD,             collapse=", "), "\n"); write.out(p.o);

  ### internal variable
}
