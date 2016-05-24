factorize <- 
function(J, D, P, omega, q) {
  vars <- unlist(lapply(P$children, FUN = function(x) x$var))
  r <- which(J == vars[omega[q]])
  J <- J[-r]
  elems <- length(J)
  if (elems > 1) {
      for (i in 1:(elems-1)) {
        j <- which(vars == J[i])
        P$children[[j]] <- probability(var = J[i], cond = union(D, J[(i+1):elems]))
      }
    }
  if (elems > 0) {
    j <- which(vars == J[elems])
    P$children[[j]] <- probability(var = J[elems], cond = D)
  }
  return(P)
}