planetFlux <-
function(T=195, dp=14.8, thetab=19.4, f=32) {
  # Function calculates the flux density from a disk-shaped planet
  #   observed in a circular Gaussian beam
  # T = planetary temperature, K
  # dp = planetary diameter, arcsec
  # thetas = beam fwhm, arcsec
  # f = frequency, GHz
  # AH 7/1/08

  c = 0.3  # units of m, GHz
  k = 1.3806e-23 # mks
  fd = 2*k*pi/(4.*log(2)) * (f/c)^2 * T * (thetab*4.848e-6)^2 *
        (1. - exp(-log(2)*(dp/thetab)^2)) * 1.e26
  return(fd)
}

