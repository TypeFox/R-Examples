#' @title List of deprecated functions
#' 
#' @description
#' These functions have been deprecated in v0.5 of \pkg{LambertW} mostly for
#'     sake of following R style guides with respect to naming of
#'     functions. This means that all deprecated functions here have an
#'     analogous function with a similar -- more style consistent -- name. See
#'     also the \code{NEWS} file.
#' 
#' Deprecated functions still work as expected, but they print out a warning
#' suggesting to use the new function (name).
#' 
#' @name deprecated-functions
#' @param ... arguments passed to deprecated functions.
#' 

NULL

#' @rdname deprecated-functions
#' @export
beta_names <- function(...) {
  
  warning("DEPRECATED: Please use get_beta_names() instead of beta_names().")
  get_beta_names(...)
}


#' @rdname deprecated-functions
#' @export
bounds_theta <- function(...) {
  warning("DEPRECATED: Please use get_theta_bounds() instead of bounds_theta().")
  get_theta_bounds(...)
}

#' @rdname deprecated-functions
#' @param z,W.z see \code{\link{deriv_W}}
#' @export
d1W_1 <- function(z, W.z = W(z, branch = -1)) {
  warning("DEPRECATED: Please usee deriv_W(z, branch = -1) instead of d1W_1(z).")
  deriv_W(z, branch = -1, W.z = W.z)
}


#' @rdname deprecated-functions
#' @export
p_1 <- function(...) {
  warning("DEPRECATED: Please use p_1m() instead of p_m1().")
  p_m1(...)
}

#' @rdname deprecated-functions
#' @export
params2theta <- function(...) {
  warning("DEPRECATED: Please use unflatten_theta() instead of unflatten_theta().")
  unflatten_theta(...)
}

#' @rdname deprecated-functions
#' @export
skewness_test <- function(...) {
  warning("DEPRECATED: Please use test_symmetry() instead of skewness_test().")
  test_symmetry(...)
}

#' @rdname deprecated-functions
#' @export
starting_theta <- function(...) {
  
  warning("DEPRECATED: Please use get_initial_theta() instead of starting_theta().")
  get_initial_theta(...)
}

#' @rdname deprecated-functions
#' @export
support <- function(...) {
  warnings("DEPRECATED: Please use 'get_support()' instead.")
  get_support(...)
}


#' @rdname deprecated-functions
#' @export
normfit <- function(...) {
  
  warning("DEPRECATED: Please use test_normality() (or 'test_norm()') ",
          "instead of normfit().")
  test_normality(...)
}


#' @rdname deprecated-functions
#' @export
theta2params <- function(...) {
  warning("DEPRECATED: use flatten_theta instead of theta2params.")
  flatten_theta(...)
}

#' @rdname deprecated-functions
#' @export
vec.norm <- function(...) {
  warning("DEPRECATED: Please use lp_norm() instead of vec.norm().")
  return(lp_norm(...))
}

#' @rdname deprecated-functions
#' @export
W_1 <- function(z) {
  # warning("DEPRECATED: Please use W(z, branch = -1) instead of W_1(z).")
  return(W(z, branch = -1))
}

#' @rdname deprecated-functions
#' @param gamma see \code{\link{W_gamma}}.
#' @export
W_gamma_1 <- function(z, gamma) {
  warning("DEPRECATED: Please use W_gamma(..., branch = -1) instead of W_gamma_1(...).")
  W_gamma(z, gamma, branch = -1)
} 


#' @rdname deprecated-functions
#' @export
H <- function(...) {
  warning("DEPRECATED: Please use xexp(...) instead of H(...).")
  return(xexp(...))
}
