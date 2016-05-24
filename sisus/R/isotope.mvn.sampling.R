isotope.mvn.sampling <-
function# Samples isotope values from the Multivariate Normal Distribution
### internal function for sisus
(n.samples.isotope.mvn
### internal variable
, isotope.mean
### internal variable
, isotope.sigma
### internal variable
)
{
  ##details<<
  ## interal function for sisus.run()

  ## removing this component to reduce package dependencies 3/5/2014
  ##library("mvtnorm");
  ##if (n.samples.isotope.mvn == 0) {
    isotope.mvn.sample = isotope.mean; # use original values if no samples drawn
    n.samples.isotope.mvn = 1;
  ##} else {
  ##  isotope.mvn.sample = rmvnorm(n.samples.isotope.mvn, mean = isotope.mean, sigma = isotope.sigma^2);
  ##}

  ISOTOPE.MVN = new.env();
  ISOTOPE.MVN$n.samples.isotope.mvn = n.samples.isotope.mvn;
  ISOTOPE.MVN$isotope.mvn.sample    = isotope.mvn.sample   ;
  return( as.list(ISOTOPE.MVN) );
  ### internal variable
}
