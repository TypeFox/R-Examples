gdbr_IBS <- 
function(gen)
{
  n = nrow(gen)
  p = ncol(gen)
  Sim = diag(1, n, n)
  # similarity distance IBS
  aux = .C("gdbr_IBS", as.integer(as.vector(t(gen))), as.integer(n), as.integer(p), as.double(as.vector(Sim)))[[4]]
  matrix(aux, nrow=n)
}
