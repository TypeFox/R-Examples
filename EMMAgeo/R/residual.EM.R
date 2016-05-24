#' Function to evaluate residual end-member loading.
#' 
#' This function calculates an optional residual end-member loading. It uses
#' the modelled end-member loadings as input and evaluates the root of 1 minus
#' the sum of all squared loadings to analyse the remaining variance, e.g.  if
#' not all (robust) EMs are included (cf. Dietze et al., 2012). Negative values
#' are set to zero.
#' 
#' 
#' @param Vqn Numeric matrix with m robust end-member loadings.
#' @return Numeric vector with residual end-member loading.
#' @author Michael Dietze, Elisabeth Dietze
#' @seealso \code{\link{EMMA}}, \code{\link{robust.EM}}
#' @references Dietze E, Hartmann K, Diekmann B, IJmker J, Lehmkuhl F, Opitz S,
#' Stauch G, Wuennemann B, Borchers A. 2012. An end-member algorithm for
#' deciphering modern detrital processes from lake sediments of Lake Donggi
#' Cona, NE Tibetan Plateau, China. Sedimentary Geology 243-244: 169-180.
#' @keywords EMMA
#' @examples
#' 
#' ## Some preparing steps to retrieve only robust end-members
#' ## load example data, i.e. here TR
#' data(rEM, envir = environment())
#' 
#' ## define mean robust end-member loadings
#' Vqn.rob <- rEM$Vqn.mean
#' ## perform residual end-member loading calculation
#' Vqn.res <- residual.EM(Vqn.rob)
#' 
#' # Visualisation of the result
#' plot(NA, xlim = c(1, 80), ylim = c(0, 1))
#' for(i in 1:4) {lines(Vqn.rob[i,])}
#' lines(Vqn.res, col = 2)
#' 
#' @export residual.EM
residual.EM <- function(
  Vqn
){
  
  ## transpose Vqn matrix
  Vqn <- t(Vqn)
  
  ## calculate squared residual end-member loading
  res.sq <- 1 - apply(Vqn^2, 1, sum)
  
  ## set negative values to zero
  res.sq[res.sq < 0] <- 0
  
  ## calculate residual end-member loading
  Vqn.res <- sqrt(res.sq)
  
  ## return result
  return(Vqn.res)
}