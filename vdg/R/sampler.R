#' Sampler Function
#' 
#' This is a wrapper for the sampling funcions of the \pkg{vdg} package. It extracts design properties from the 
#' design passed to it to take appropriate samples.
#' 
#' @param n number of points to sample
#' @param design  design for which the sample is required (either a matrix or data frame)
#' @param type type of design region/sampling method. One of "spherical", "cuboidal", 
#' "lhs", "mlhs", "slhs", "rslhs" or "custom". Option "custom" requires \code{custom.fun} to be
#' non-\code{NULL}.
#' @param at logical; should sampling be done on the surface of hyperspheres or hypercubes? 
#' Not used for LHS methods.
#' @param custom.fun A custom sampling function, used in conjunction with \code{type = "custom"}. The
#' first and second arguments must be the sample size and dimension respectively.
#' @param \dots other arguments passed to the underlying sampling functions.
#' @return Matrix with samples as rows, with S3 class \code{smpl}
#' @seealso \code{\link{runif_sphere}}, \code{\link{runif_cube}}, \code{\link{LHS}}, 
#' \code{\link{MLHS}}, \code{\link{SLHS}}, \code{\link{RSLHS}}
#' @author Pieter C. Schoonees
#' @export
#' @examples
#' ## Default spherical design region
#' set.seed(1896)
#' samp1 <- sampler(n = 100, design = expand.grid(x = -1:1, y = -1:1))
#' plot(samp1)
#' 
#' ## Supplying a custom sampling function based on lhs::improvedLHS()
#' library("lhs")
#' sfun <- function(n, k, dup = 1) 2 * improvedLHS(n, k, dup = dup) - 1
#' samp2 <- sampler(n = 100, design = expand.grid(x = -1:1, y = -1:1),
#'                  type = "custom", custom.fun = sfun)
#' plot(samp2)
sampler <- function(n, design, type = c("spherical", "cuboidal", "lhs", "mlhs", "slhs", "rslhs", "custom"), 
                    at = FALSE, custom.fun = NULL, ...){
  m <- ncol(design) 
  type <- tolower(type)
  type <- match.arg(arg = type)
  if (type == "custom" && !is.function(custom.fun)) 
    stop("A custom sampling function must be supplied as argument 'custom.fun'.")
  samp <- switch(type, spherical = runif_sphere(n = n, m = m, at = at, ...), 
                 cuboidal = runif_cube(n = n, m = m, at = at, ...), 
                 lhs = LHS(n = n, m = m, ...), 
                 mlhs = MLHS(n = n, m = m, ...),
                 slhs = SLHS(n = n, m = m, ...), 
                 rslhs = RSLHS(n = n, m = m, ...),
                 custom = custom.fun(n, m, ...))
  colnames(samp) <- colnames(design)
  class(samp) <- c("smpl", "matrix")
  samp
}