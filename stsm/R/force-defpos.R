
force.defpos <- function(m, tol = 0.001, debug = FALSE)
#
# matrix 'm' is forced to be definite positive by
# changing the eigenvalues
# Pollock:99 pp.342, Goldfeld-etal
# see also other alternatives in Nocedal-Wright 2006 Chapter 3
#
{
  #is.Symmetric(IM)
  eg <- eigen(m, only.values = TRUE)$values #symmetric = TRUE
  min.eg <- min(eg)
  if (min.eg <= 0)
    m <- m + (tol - min.eg) * diag(nrow(m))

  if (debug)
  {
    eg <- eigen(m, only.values = TRUE)$values
    if (min(eg) <= 0) {
      warning("the matrix could not be rendered into a positive definite matrix.")
      #print(eg)
    }
  }

  m
}
