emul.minus.lik <-
function(parvec, Y.mat, X.mat, t.vec, Theta.mat, n.par, p.par, fix.betas,
                     limits.lower=NULL, limits.upper=NULL, beta.vec=NULL) { 
  llik <- emul.lik(parvec, Y.mat, X.mat, t.vec,
             Theta.mat, n.par, p.par, fix.betas, limits.lower, limits.upper, beta.vec)
  mllik <- (-1)*llik 
  mllik 
}
