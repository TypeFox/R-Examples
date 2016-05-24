prior_unif_dirichlet <-
function(p,r,Neq = 0.5*r^2){
  a = r*diag(numeric(r)+1)
  b = matrix(1,r,r)
  return((Neq/(r^2))*array(c(a,rep(b,p)),c(r,r,p,p)))
}
