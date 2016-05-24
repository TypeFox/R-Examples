
## non-exported function for low-level testing
.gaussian_near_field <- function(xyz, wavelength=500, alpha = 15*pi/180, psi=0,
                                 w0=1e4, ni=1.5, no=1.0, nl=no, d=0,
                                 cutoff = min(1, 3*wavelength/(ni*pi*w0)), # integrand ~ exp(-3^2)
                                 maxEval = 3000, tol=1e-04, field=FALSE){

  k0 <- 2*pi/wavelength
  ki <- k0*ni
  res <- cubature::adaptIntegrate(cpp_integrand_gb_layer,
                                  lowerLimit=c(0, 0), # rho in [0,1], angle in [0,2*pi]
                                  upperLimit=c(cutoff, 2*pi),
                                  fDim = 6, tol = tol,
                                  maxEval = maxEval,
                                  r2 = xyz, ki=ki, psi= psi, alpha=alpha,
                                  w0=w0, ni=ni, no=no, nl=nl, d=d)$integral

  E <- complex(real = res[c(1,3,5)], imaginary=res[c(2,4,6)])

  if(field) return(E)

  Re(crossprod(E, Conj(E)))

}

##' Electric field from the transmission of a gaussian beam at a planar interface
##'
##' Integration is performed over a spectrum of incident plane waves
##' @title gaussian_near_field_layer
##' @param xyz position
##' @param wavelength wavelength
##' @param alpha beam incident angle
##' @param psi beam polarisation angle
##' @param w0 beam waist radius
##' @param epsilon vector of permittivities
##' @param thickness thickness corresponding to each medium
##' @param maxEval passed to adaptIntegrate
##' @param reqAbsError passed to cubature
##' @param tol passed to adaptIntegrate
##' @param progress logical: display progress bar
##' @param field logical: return the electric field (complex vector), or modulus squared
##' @return data.frame electric field at the x, y, z position
##' @export
##' @family gaussian_beam
##' @author Baptiste Auguie
gaussian_near_field_layer <- function(xyz, wavelength=500, alpha = 15*pi/180, psi=0,
                           w0=1e4, epsilon = c(1.5^2, epsAg(wavelength)$epsilon, 1.0^2),
                           thickness = c(0, 50, 0),
                           maxEval = 3000, reqAbsError=0.0, tol=1e-04,
                           progress = FALSE, field=FALSE){

  k0 <- 2*pi/wavelength
  E <- cpp_field_gb_layer(xyz, k0, psi, alpha, w0, epsilon, thickness,
                          as.integer(maxEval), as.double(reqAbsError), tol, progress)

  if(field) return(E)

  Re(colSums(E*Conj(E)))

}



##' Electric field of a gaussian beam close to a planar interface
##'
##' Integration is performed over a spectrum of incident plane waves using integrand_gb2
##' @title gaussian_near_field_ml
##' @param xyz position matrix
##' @param wavelength wavelength
##' @param alpha beam incident angle
##' @param psi beam polarisation angle
##' @param w0 beam waist radius
##' @param epsilon vector of permittivities
##' @param thickness thickness corresponding to each medium
##' @param maxEval passed to cubature
##' @param reqAbsError passed to cubature
##' @param tol passed to cubature
##' @param progress logical display progress bar
##' @param field logical: return the electric field (complex vector), or modulus squared
##' @return data.frame electric field at the x, y, z position
##' @export
##' @family gaussian_beam
##' @author Baptiste Auguie
gaussian_near_field_ml <- function(xyz, wavelength=632.8, alpha = 15*pi/180, psi=0,
                                 w0=1e4, epsilon = c(1.5^2, epsAg(wavelength)$epsilon, 1.0^2, 1.0^2),
                                 thickness = c(0, 50, 10, 0),
                                 maxEval = 3000, reqAbsError = 0.0, tol=1e-04, progress=FALSE, field=FALSE){

  k0 <- 2*pi/wavelength
  E <- cpp_field_gb_ml(xyz, k0, psi, alpha, w0, epsilon, thickness,
                         as.integer(maxEval), as.double(reqAbsError), tol, progress)

  if(field) return(E)

  Re(colSums(E*Conj(E)))
}
