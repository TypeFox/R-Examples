# utility function for combining two fgh lists by adding their
# corresponding elements; useful for adding log-prior and log-likelihood
regfac.merge <- function(fgh1, fgh2, fgh=2) {
  if (fgh==0) return (fgh1+fgh2)
  return (list(f=fgh1$f+fgh2$f, g=fgh1$g+fgh2$g, h=fgh1$h+fgh2$h))
}

