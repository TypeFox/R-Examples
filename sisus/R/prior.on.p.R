prior.on.p <-
function# Prior on the vector of source proportions, p, is defined by the given vector alpha
### internal function for sisus
(priors.sources
### internal variable
, priors.precision
### internal variable
, n.sources
### internal variable
)
{
  ##details<<
  ## interal function for sisus.run()

  alpha.p.dirichlet = n.sources * priors.sources;
  dirichlet.moments(alpha.p.dirichlet, priors.precision, priors.sources);
  return(alpha.p.dirichlet);
  ### internal variable
}
