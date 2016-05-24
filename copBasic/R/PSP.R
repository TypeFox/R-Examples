"PSP" <-
function(u,v,...) {
  if(length(u) > 1 & length(v) > 1 & length(u) != length(v)) {
    warning("length u = ",length(u), " and length v = ",length(v))
    warning("longer object length is not a multiple of shorter object length, ",
            "no recycling in PSP()")
    return(NA)
  }
  # The extra hassle of vectorization made here is to handle situations
  # in which nested integrals are used where uneven vectors can be passed
  # Situation came to life years after initial release in testing
  # rhoCOP(W) and rhoCOP(M) for which Nelsen (2006) provides analytical solutions
  if(length(u) == 1) {
     u <- rep(u, length(v))
  } else if(length(v) == 1) {
     v <- rep(v, length(u))
  }
  return(sapply(1:length(u), function(i) {
                    exp(log(u[i]) + log(v[i]) - log(u[i] + v[i] - (u[i]*v[i]))) }))
}
