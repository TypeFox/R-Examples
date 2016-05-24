######################################################################
# Sample size based on Wald score confidence intervals (the usual thing)
#
# Params:
#   p0   - hypothesized upper bound (if below 0.5, if above 0.5 then
#          lower bound) on the parameter p.
#   d    - half width of the confidence interval
#  alpha - level of the confidence interval
######################################################################

ciss.wald <- function(p0, d, alpha) {
  if (d >= 1/2) {
    stop("Error: d has to be smaller than 1/2.")
  }
  
  #Compute values used in Eq (5) of Piegorsch (2004)
  z2 <- qnorm(alpha/2)^2
  return(ceiling( z2*p0*(1-p0)/d^2))
}


######################################################################
# Sample size based on Wilson score confidence intervals
#
# Params:
#   p0   - hypothesized upper bound (if below 0.5, if above 0.5 then
#          lower bound) on the parameter p.
#   d    - half width of the confidence interval
#  alpha - level of the confidence interval
######################################################################

ciss.wilson <- function(p0, d, alpha) {
  if (d >= 1/2) {
    stop("Error: d has to be smaller than 1/2.")
  }
  
  #Compute values used in Eq (7) of Piegorsch (2004)
  z2 <- qnorm(alpha/2)^2
  d2 <- d^2
  
  nS <- z2/(2*d2) * (p0 * (1-p0) - 2*d2 + sqrt( p0^2*(1-p0)^2 + 4*d2 * (p0-1/2)^2 ))

  return(ceiling(nS))
}

######################################################################
# Sample size based on the Agresti-Coull confidence intervals
#
# Params:
#   p0   - hypothesized upper bound (if below 0.5, if above 0.5 then
#          lower bound) on the parameter p.
#   d    - half width of the confidence interval
#  alpha - level of the confidence interval
######################################################################

ciss.agresticoull <- function(p0, d, alpha) {
  if (d >= 1/2) {
    stop("Error: d has to be smaller than 1/2.")
  }
  #Compute according to equation (6) in Piegorsch (2004)
  z2 <-  qnorm(alpha/2)^2
  n <- z2*p0*(1-p0)/d^2 - z2
  return(ceiling(n))
}

