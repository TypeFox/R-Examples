probability <-
function(var = character(), cond = character(), sumset = character(), do = character(), recursive = FALSE, children = list(), fraction = FALSE, divisor = list(), star = FALSE, domain = character()) {
  p <- list(var = var, cond = cond, sumset = sumset, do = do, recursive = recursive, children = children, fraction = fraction, divisor = divisor, star = star, domain = domain)
  class(p) = "probability"
  return(p)
}
