##' C++ calculation of the array factor
##'
##' Brute-force numerical evaluation of the truncated 2D sum of dipole fields
##' @title array factor
##' @param wavelength wavelength in nm
##' @param N half number of dipoles along one side
##' @param pitch pitch in nm
##' @return S
##' @export
##' @family user_level array
##' @author baptiste Auguie
##'  @examples
##' S <- array_factor(seq(400, 600),  10,  500)
##' str(S)
array_factor <- function(wavelength, N, pitch){

  k <- 2*pi/wavelength
  rj <- expand.grid(-N:N, -N:N) * pitch
  rj <- rj[-((dim(rj)[1]-1)/2 + 1),1:2] # remove rj=(0,0)
  S <- array$array_factor(k, as.matrix(rj))

  invisible(data.frame(wavelength=wavelength, S=S))
  
}
