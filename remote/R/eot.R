if ( !isGeneric('eot') ) {
  setGeneric('eot', function(x, ...)
    standardGeneric('eot'))
}

#' EOT analysis of a predictor and (optionally) a response RasterStack
#' 
#' @description
#' Calculate a given number of EOT modes either internally or between 
#' RasterStacks.
#' 
#' @param x a RasterStack used as predictor
#' @param y a RasterStack used as response. If \code{y} is \code{NULL},
#' \code{x} is used as \code{y}
#' @param n the number of EOT modes to calculate
#' @param standardised logical. If \code{FALSE} the calculated r-squared values 
#' will be multiplied by the variance
#' @param write.out logical. If \code{TRUE} results will be written to disk 
#' using \code{path.out}
#' @param path.out the file path for writing results if \code{write.out} is \code{TRUE}.
#' Defaults to current working directory
#' @param prefix optional prefix to be used for naming of results if 
#' \code{write.out} is \code{TRUE}
#' @param reduce.both logical. If \code{TRUE} both \code{x} and \code{y} 
#' are reduced after each iteration. If \code{FALSE} only \code{y} is reduced
#' @param type the type of the link function. Defaults to \code{'rsq'} as in original
#' proposed method from \cite{van den Dool 2000}. If set to \code{'ioa'} index of agreement is
#' used instead
#' @param verbose logical. If \code{TRUE} some details about the 
#' calculation process will be output to the console
#' @param ... not used at the moment
#' 
#' @details 
#' For a detailed description of the EOT algorithm and the mathematics behind it,
#' see the References section. In brief, the algorithm works as follows: 
#' First, the temporal profiles of each pixel \emph{xp} of the predictor domain 
#' are regressed against the profiles of all pixels \emph{xr} in the 
#' response domain. 
#' The calculated coefficients of determination are summed up and the pixel 
#' with the highest sum is identified as the 'base point' of the first/leading mode. 
#' The temporal profile at this base point is the first/leading EOT. 
#' Then, the residuals from the regression are taken to be the basis 
#' for the calculation of the next EOT, thus ensuring orthogonality 
#' of the identified teleconnections. This procedure is repeated until 
#' a predefined amount of \emph{n} EOTs is calculated. In general, 
#' \pkg{remote} implements a 'brute force' spatial data mining approach to 
#' identify locations of enhanced potential to explain spatio-temporal 
#' variability within the same or another geographic field.
#' 
#' @return 
#' if n = 1 an \emph{EotMode}, if n > 1 an \emph{EotStack} of \code{n} 
#' \emph{EotMode}s. Each \emph{EotMode} has the following components:
#' \itemize{
#' \item \emph{mode} - the number of the identified mode (1 - n)
#' \item \emph{eot} - the EOT (time series) at the identified base point. 
#' Note, this is a simple numeric vector, not of class \code{ts}
#' \item \emph{coords_bp} - the coordinates of the identified base point
#' \item \emph{cell_bp} - the cell number of the indeified base point
#' \item \emph{cum_exp_var} - the (cumulative) explained variance of the considered EOT
#' \item \emph{r_predictor} - the \emph{RasterLayer} of the correlation coefficients 
#' between the base point and each pixel of the predictor domain
#' \item \emph{rsq_predictor} - as above but for the coefficient of determination
#' \item \emph{rsq_sums_predictor} - as above but for the sums of coefficient of determination
#' \item \emph{int_predictor} - the \emph{RasterLayer} of the intercept of the 
#' regression equation for each pixel of the predictor domain
#' \item \emph{slp_predictor} - same as above but for the slope of the 
#' regression equation for each pixel of the predictor domain
#' \item \emph{p_predictor} - the \emph{RasterLayer} of the significance (p-value) 
#' of the the regression equation for each pixel of the predictor domain
#' \item \emph{resid_predictor} - the \emph{RasterBrick} of the reduced data 
#' for the predictor domain
#' }
#' 
#' Apart from \emph{rsq_sums_predictor}, all \emph{*_predictor} fields are 
#' also returned for the \emph{*_response} domain, 
#' even if predictor and response domain are equal. This is due to that fact, 
#' that if not both fields are reduced after the first EOT is found, 
#' these \emph{RasterLayers} will differ.
#' 
#' 
#' @references 
#' \bold{Empirical Orthogonal Teleconnections}\cr
#' H. M. van den Dool, S. Saha, A. Johansson (2000)\cr
#' Journal of Climate, Volume 13, Issue 8, pp. 1421-1435\cr
#' \url{http://journals.ametsoc.org/doi/abs/10.1175/1520-0442%282000%29013%3C1421%3AEOT%3E2.0.CO%3B2
#' }
#'  
#' \bold{Empirical methods in short-term climate prediction}\cr
#' H. M. van den Dool (2007)\cr
#' Oxford University Press, Oxford, New York\cr
#' \url{http://www.oup.com/uk/catalogue/?ci=9780199202782}
#' 
#' @examples
#' ### EXAMPLE I
#' ### a single field
#' data(vdendool)
#' 
#' ## claculate 2 leading modes
#' nh_modes <- eot(x = vdendool, y = NULL, n = 2, 
#'                 reduce.both = FALSE, standardised = FALSE, 
#'                 verbose = TRUE)
#' 
#' plot(nh_modes, y = 1, show.bp = TRUE)
#' plot(nh_modes, y = 2, show.bp = TRUE)
#' 
#' @export
#' @name eot
#' @rdname eot
#' @aliases eot,RasterStack-method

# set methods -------------------------------------------------------------

setMethod('eot', signature(x = 'RasterStack'), 
          function(x, 
                   y = NULL, 
                   n = 1, 
                   standardised = TRUE, 
                   write.out = FALSE,
                   path.out = ".", 
                   prefix = "remote",
                   reduce.both = FALSE, 
                   type = c("rsq", "ioa"),
                   verbose = TRUE,
                   ...) {
            
            # Duplicate predictor set in case predictor and response are identical
            if (is.null(y)) {
              y <- x  
            }
            
            orig.var <- calcVar(y, standardised = standardised)
            
            ### EOT
            
            # Loop through number of desired EOTs
            for (z in 1:n) {
              
              # Use initial response data set in case of first iteration
              if (z == 1) {
                
                x.eot <- EotCycle(x = x, 
                                  y = y,
                                  n = z, 
                                  type = type,
                                  standardised = standardised, 
                                  orig.var = orig.var,
                                  write.out = write.out,
                                  path.out = path.out, 
                                  verbose = verbose,
                                  prefix = prefix)
                
                names(x.eot) <- paste("mode_", sprintf("%02.f", z), 
                                      sep = "")
                
                # Use last entry of slot 'residuals' otherwise  
              } else if (z > 1) {
                tmp.x.eot <- EotCycle(
                  x = if (!reduce.both) {
                    x
                  } else {
                    if (z == 2) {
                      x.eot@resid_predictor
                    } else {
                      x.eot[[z-1]]@resid_predictor
                    }
                  }, 
                  y = if (z == 2) {
                    x.eot@resid_response 
                  } else {
                    x.eot[[z-1]]@resid_response
                  }, 
                  y.eq.x = y.eq.x,
                  n = z, 
                  type = type,
                  standardised = standardised, 
                  orig.var = orig.var,
                  write.out = write.out,
                  path.out = path.out,  
                  verbose = verbose,
                  prefix = prefix)
                
                if (z == 2) {
                  x.eot <- list(x.eot, tmp.x.eot)
                  names(x.eot) <- c(paste("mode_", sprintf("%02.f", 1), 
                                          sep = ""), 
                                    paste("mode", sprintf("%02.f", z), 
                                          sep = "_"))
                } else {
                  tmp.names <- names(x.eot)
                  x.eot <- append(x.eot, list(tmp.x.eot))
                  names(x.eot) <- c(tmp.names, 
                                    paste("mode", sprintf("%02.f", z), 
                                          sep = "_"))
                }
              }
            }
            
            if (length(x.eot) == 1) {
              out <- x.eot
            } else {
              out <- new('EotStack', modes = x.eot, names = names(x.eot))
            }
            return(out)
          }
)

#' @describeIn eot

setMethod('eot', signature(x = 'RasterBrick'), 
          function(x, 
                   y = NULL, 
                   n = 1, 
                   standardised = TRUE, 
                   write.out = FALSE,
                   path.out = ".", 
                   prefix = "remote",
                   reduce.both = FALSE, 
                   type = c("rsq", "ioa"),
                   verbose = TRUE,
                   ...) {
            
            # Duplicate predictor set in case predictor and response are identical
            if (is.null(y)) {
              y <- x  
              y.eq.x <- TRUE
            } else {
              y.eq.x <- FALSE
            }
            
            orig.var <- calcVar(y, standardised = standardised)
            
            ### EOT
            
            # Loop through number of desired EOTs
            for (z in seq(n)) {
              
              # Use initial response data set in case of first iteration
              if (z == 1) {
                
                x.eot <- EotCycle(x = x, 
                                  y = y,
                                  y.eq.x = y.eq.x,
                                  n = z, 
                                  type = type,
                                  standardised = standardised, 
                                  orig.var = orig.var,
                                  write.out = write.out,
                                  path.out = path.out, 
                                  verbose = verbose,
                                  prefix = prefix)
                
                names(x.eot) <- paste("mode_", sprintf("%02.f", z), 
                                      sep = "")
                
                # Use last entry of slot 'residuals' otherwise  
              } else if (z > 1) {
                tmp.x.eot <- EotCycle(
                  x = if (!reduce.both) {
                    x
                  } else {
                    if (z == 2) {
                      x.eot@resid_predictor
                    } else {
                      x.eot[[z-1]]@resid_predictor
                    }
                  }, 
                  y = if (z == 2) {
                    x.eot@resid_response 
                  } else {
                    x.eot[[z-1]]@resid_response
                  }, 
                  y.eq.x = y.eq.x,
                  n = z, 
                  type = type,
                  standardised = standardised, 
                  orig.var = orig.var,
                  write.out = write.out,
                  path.out = path.out,  
                  verbose = verbose,
                  prefix = prefix)
                
                if (z == 2) {
                  x.eot <- list(x.eot, tmp.x.eot)
                  names(x.eot) <- c(paste("mode_", sprintf("%02.f", 1), 
                                          sep = ""), 
                                    paste("mode", sprintf("%02.f", z),
                                          sep = "_"))
                } else {
                  tmp.names <- names(x.eot)
                  x.eot <- append(x.eot, list(tmp.x.eot))
                  names(x.eot) <- c(tmp.names, 
                                    paste("mode", sprintf("%02.f", z), 
                                          sep = "_"))
                }
              }
            }
            
            if (length(x.eot) == 1) {
              out <- x.eot
            } else {
              out <- new('EotStack', modes = x.eot, names = names(x.eot))
            }
            return(out)
          }
)

