wkde <- function(x, u=x, w=rep(1, length(x)), bw=bw.nrd0(as.vector(x)), sym=FALSE) {
  if (sym) {
    return((wkde(x, u, w, bw) + wkde(x, -u, w, bw))/2)
  }
  Km <- exp(outer(x/bw, u/bw, function(a,b) -(a-b)^2/2))
	normw <- matrix(w/sum(w), nrow=1)
  as.vector(normw %*% Km) / (bw * sqrt(2*pi))
}

