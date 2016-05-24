if (!isGeneric('nXplain')) {
  setGeneric('nXplain', function(x, ...)
    standardGeneric('nXplain')) 
}


#' Number of EOTs needed for variance explanation
#' 
#' @description 
#' The function identifies the number of modes needed to explain a certain amount of
#' variance within the response field.
#' 
#' @param x an \emph{EotStack}
#' @param var the minimum amount of variance to be explained by the modes
#' 
#' @note This is a post-hoc function. It needs an \emph{EotStack} 
#' created as returned by \code{\link{eot}}. Depending on the potency
#' of the identified EOTs, it may be necessary to compute a high number of 
#' modes in order to be able to explain a large enough part of the variance.
#' 
#' @return an integer denoting the number of EOTs needed to explain \code{var}
#' 
#' @examples
#' data(vdendool)
#' 
#' nh_modes <- eot(x = vdendool, y = NULL, n = 3, 
#'                 reduce.both = FALSE, standardised = FALSE, 
#'                 verbose = TRUE)
#'              
#' ### How many modes are needed to explain 25% of variance?              
#' nXplain(nh_modes, 0.25)
#' 
#' @export 
#' @name nXplain
#' @rdname nXplain
#' @aliases nXplain,EotStack-method

setMethod('nXplain', signature(x = 'EotStack'),
          function(x, var = 0.9) {
            expl.var <- sapply(seq(nmodes(x)), function(i) {
              x[[i]]@cum_exp_var
            })
            
            n <- min(which(var - expl.var <= 0), na.rm = TRUE)
            
            if (!is.finite(n)) {
              stop("explained variance of EotStack is lower than: ", var, 
                   "\n",
                   "maximum explained variance of this EotStack is: ", 
                   x[[nmodes(x)]]@cum_exp_var)
            }
            
            return(n)
          }
)