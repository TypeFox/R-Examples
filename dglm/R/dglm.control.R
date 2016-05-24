dglm.control <- function(epsilon = 1e-007, maxit = 50, trace = FALSE, ...) {
  #  Control iteration for dglm
  #  As for glm.control but with different defaults
  #  GKS  9 Jan 98
  #
  if (epsilon <= 0) {
    warning("the value of epsilon supplied is zero or negative;
            the default value of 1e-7 was used instead")
    epsilon <- 1e-007}
  if (maxit < 1) {
    warning("the value of maxit supplied is zero or negative;
            the default value of 50 was used instead")
    maxit <- 50}
  list(epsilon = epsilon, maxit = maxit, trace = as.logical(trace)[1])
}
