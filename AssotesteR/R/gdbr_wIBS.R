gdbr_wIBS <- 
function(gen, weights)
{
  n = nrow(gen)
  p = ncol(gen)
  Sim = diag(1, n, n)
  aux = .C("gdbr_wIBS", as.integer(as.vector(t(gen))), as.integer(n), as.integer(p), 
           as.double(weights), as.double(as.vector(Sim)))[[5]]
  matrix(aux, nrow=n)
}
