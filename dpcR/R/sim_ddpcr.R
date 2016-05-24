#' Simulate Droplet Digital PCR
#' 
#' A function that simulates results of a droplet digital PCR.
#' 
#' The function contains two
#' implementations of the array digital PCR simulation. First one was described
#' in Dube at. al (2008). This method is based on random distributing \eqn{m
#' \times times}{m * times} molecules between \eqn{n \times times}{n * times}
#' chambers.  After this step, the required number of plates is created by the
#' random sampling of chambers without replacement. The above method is used,
#' when the \code{dube} argument has value \code{TRUE}.
#' 
#' The higher the value of the argument \code{times}, the simulation result is
#' closer to theoretical calculations.
#' 
#' @param m the total number of template molecules used in the expertiment. Must be
#' a positive integer.
#' @param n the number of droplets per experiment. Must be a positive integer.
#' @param times number of repetitions (see Details).
#' @param n_exp the number of experiments that are simulated by the function.
#' Cannot have higher value than the \code{times} argument.
#' @param dube if \code{TRUE}, the function is strict implementation of digital
#' PCR simulation (as in Dube et al., 2008). If \code{FALSE}, the function
#' calculates only approximation of Dube's experiment. See Details and
#' References.
#' @param pos_sums if \code{TRUE}, function returns only the total number of
#' positive (containing at least one molecule) chamber per panel. If
#' \code{FALSE}, the functions returns a vector of length equal to the number
#' of chambers. Each element of the vector represents the number of template
#' molecules in a given chamber.
#' @param fluo if \code{NULL}, the function calculates number of molecules per
#' well or total number of positive droplets. If list of two, the first
#' argument defines smoothness of the fluorescence curve and second space
#' between two consecutive measured droplets. Space must be a vector containing
#' positive integers of the length \code{n} or 1.
#' @return If the \code{pos_sums} argument has value \code{FALSE}, the function
#' returns matrix with \eqn{n} rows and \eqn{n_panels} columns. Each column
#' represents one plate. The type of such simulation would be \code{"nm"}. If the
#' \code{pos_sums} argument has value \code{TRUE}, the function return matrix
#' with one row and \eqn{n_panels} columns. Each column contains the total
#' number of positive chambers in each plate and type of simulation would be
#' set as \code{"tnp"}.
#' 
#' In each case the value is an object of the \code{\linkS4class{ddpcr}} class.
#' @note Although Dube's simulation of digital PCR was developed for array
#' digital PCR, it's also viable for simulating droplet-based methods.
#' @author Michal Burdukiewicz, Stefan Roediger.
#' @seealso \code{\link{sim_adpcr}}.
#' @keywords datagen
#' @examples
#' 
#' #simulate fluorescence data
#' tmp_VIC <- sim_ddpcr(m = 7, n = 20, times = 5, fluo = list(0.1, 0))
#' tmp_FAM <- sim_ddpcr(m = 15, n = 20, times = 5, fluo = list(0.1, 0))
#' par(mfrow = c(2,1))
#' plot(tmp_VIC, col = "green", type = "l")
#' plot(tmp_FAM, col = "blue", type = "l")
#' summary(tmp_FAM)
#' 
#' summary(sim_ddpcr(m = 7, n = 20, times = 5, n_exp = 5))
#' 
#' @export sim_ddpcr
sim_ddpcr <- function(m, n, times, n_exp = 1, dube = FALSE, pos_sums = FALSE, 
                      fluo = NULL) {
  if (!is.null(fluo))
    if (pos_sums)
      stop("During fluorescence simulation 'pos_sums' must be TRUE", call. = TRUE, 
           domain = NA)
  n <- num2int(n)
  res <- sim_dpcr(m, n, times, dube, pos_sums, n_exp)
  if (!is.null(fluo)) {
    res <- apply(res, 2, function(x) sim_ddpcr_fluo(x, n, fluo[[1]], fluo[[2]]))
  }
  #simplify
  type = ifelse(pos_sums, "tnp", "nm")
  if (!is.null(fluo))
    type <- "fluo"
  create_ddpcr(res, n = rep(n, n_exp), threshold = 0.5, type = type)
}


sim_ddpcr_fluo <- function(res, n, resolution, space) {
  if (length(space) > 1) {
    a <- lapply(space, function(x) c(sin(seq(0, pi, resolution)), rep.int(0, x)))
  } else {
    a <- lapply(1L:n, function(x) c(sin(seq(0, pi, resolution)), rep.int(0, space)))
  }
  result <- matrix(unlist(lapply(1L:n, function(x)
    a[[x]] * res[x])), ncol = 1)
  result
}
