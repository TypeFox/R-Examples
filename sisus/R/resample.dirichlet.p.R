resample.dirichlet.p <-
function# Resamples weighted samples from the polytope from the Dirichlet distribution on the p vector
### internal function for sisus
(M.re
### internal variable
, sam.re
### internal variable
, alpha.p.dirichlet
### internal variable
, n.vertices
### internal variable
)
{
  ##details<<
  ## interal function for sisus.run()

  # attached via DESCRIPTION # library("gtools");
  # Dirichlet resampling of p-vectors
  # Weights
  weights.p.dirichlet.temp = ddirichlet(sam.re[(n.vertices+1):M.re,], alpha.p.dirichlet);             # evaluate dirichlet probabilities
  weights.p.dirichlet = weights.p.dirichlet.temp/sum(weights.p.dirichlet.temp);                       # make weights probabilities
  # resample
  resample.ind = sample(seq(n.vertices+1, M.re), M.re, replace = TRUE, prob = weights.p.dirichlet);   # draw sample indicies
  resam = rbind(sam.re[1:n.vertices,], sam.re[resample.ind,]);                                        # update sampling with resample, putting vertices at beginning of sample
  if (M.re==1) {resam = t(as.matrix(resam));}                                                         # if single solution, still make a matrix
  return(resam);                                                                                      # return resample
  ### internal variable
}
