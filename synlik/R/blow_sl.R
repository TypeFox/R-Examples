#' Blowfly model
#' 
#' @description \code{synlik} object containing the blowfly model proposed by Wood (2010).
#'              The main components are the simulator \link{blowSimul} and the statistics
#'              \code{blowStats}, described in the same reference.
#'  
#'
#' @name blow_sl
#' @rdname blow_sl
#' @aliases blowStats blow_smcmc
#' @seealso \link{blowSimul}
#' @export
#' @references Simon N Wood. Statistical inference for noisy nonlinear ecological dynamic systems. Nature, 466(7310):1102--1104, 2010.
#' @author Simon Wood and Matteo Fasiolo <matteo.fasiolo@@gmail.com>.
#' @examples
#' data(blow_sl)
#' 
#' plot(blow_sl)
#' simulate(blow_sl, stats = TRUE)
#' 
#' slik(blow_sl, 
#'      param  = log( c( "delta" = 0.16, "P" = 6.5, "N0" = 400, 
#'                       "var.p" = 0.1, "tau" = 14, "var.d" = 0.1)  ),
#'      nsim    = 1e3)
#' 
#' # Using Nicholson's data
#' data(bf1)
#' 
#' plot(blow_sl)
#' 
#' blow_sl@@data <- bf1$pop
#' blow_sl@@extraArgs$obsData <- bf1$pop #Important: blow_sl@@blowStats uses the observed data
#' 
#' slik(blow_sl, 
#'      param  = log( c( "delta" = 0.16, "P" = 6.5, "N0" = 400, 
#'                       "var.p" = 0.1, "tau" = 14, "var.d" = 0.1)  ),
#'      nsim    = 1e3)
NULL