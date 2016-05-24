##' Drude model for the dielectric function of good (governed by free electrons) metals
##'
##' a bacground contribution eps_inf is assumed for the core electrons
##' @title drude
##' @param wavelength wavelength in nm
##' @param p vector of 3 parameters
##' @param omega angular frequency in rad/s
##' @param omega_p plasma frequency in rad/s
##' @param gamma_p damping constant, in rad/s
##' @param epsilon_inf background dielectric function
##' @param ... not used
##' @return a data.frame with wavelength in nm and complex dielectric function
##' @author Baptiste Auguie
##' @export
drude <- function(wavelength=633, p=c(1e16, 1e14, 1), 
                  omega=2*pi*1e9 * 2.99792458e8 / wavelength, 
                  omega_p = p[1], gamma_p = p[2], epsilon_inf = p[3], ...){
  
  if(missing(wavelength))
    wavelength <- 2*pi*1e9 * 2.99792458e8 / omega
  
  data.frame(wavelength=wavelength,
             epsilon = epsilon_inf -  omega_p^2 / (omega^2 +  (0+1i) * omega * gamma_p))
  
}

##' Objective function for the Drude model
##'
##' Used to fit a Drude model to a material
##' @title fit_drude
##' @param p parameters vector (3)
##' @param material data.frame with wavelength in nm and complex epsilon
##' @param ... passed to drude
##' @return sum of squares
##' @author Baptiste Auguie
##' @export
fit_drude <- function(p, material, ...){
  
  res <- drude(material[['wavelength']], p, ...)
  
  sos.real <- sum((Re(res[['epsilon']]) - Re(material[['epsilon']]))^2) / 
    sum(Re(material[['epsilon']])^2)
  sos.imag <- sum((Im(res[['epsilon']]) - Im(material[['epsilon']]))^2) / 
    sum(Im(material[['epsilon']])^2)
  
  sos.real + sos.imag
}
