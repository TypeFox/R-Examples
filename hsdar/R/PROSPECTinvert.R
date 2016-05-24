.inversionFUN <- function(P, reflectance_spectra, transmittance_spectra, sam)
{
  pros_spectra <- spectra(PROSPECT(N = P[1]/10, Cab = P[2], Car = P[3], Cbrown = P[4]/100,
                                   Cw = P[5]/100, Cm = P[6]/100, transmittance = FALSE,
                                   parameterList = NULL))[1,]
  pros_spectra[!is.finite(pros_spectra)] <- 0                                   
  if (is.null(transmittance_spectra))
  {
    
    if (sam)
    {
      specang <- 0
      chi2 <- .Fortran("sam",
                       nspec=as.integer(1),
                       nref=as.integer(1),
                       nbands=as.integer(2101),
                       spec=as.double(pros_spectra),
                       specref=as.double(reflectance_spectra),
                       specang=as.double(specang),
                       PACKAGE="hsdar"
                      )$specang
    } else {
      chi2 <- sqrt(mean((reflectance_spectra-pros_spectra)^2))
    }
  } else {
    pros_spectra_t <- spectra(PROSPECT(N = P[1]/10, Cab = P[2], Car = P[3], Cbrown = P[4]/100,
                                       Cw = P[5]/100, Cm = P[6]/100, transmittance = FALSE,
                                       parameterList = NULL))[1,]
    chi2 <- sqrt(sum((pros_spectra-reflectance_spectra)^2+(pros_spectra_t-transmittance_spectra)^2))
  }
  return(chi2)
}
  
PROSPECTinvert <- function(x, P0 = NULL, transmittance_spectra = NULL, sam = FALSE, ...)
{
  if (!requireNamespace("pracma", quietly = TRUE))
    stop("Library 'pracma' is required to invert PROSPECT")
  
  if (is.null(P0))
  {
    P0 <- (c(10, 0, 0, 0, 0.005, 0.2) + c(30, 100, 30, 100, 4, 1.8))/2
  } else {
    P0 <- c(P0[1]*10, P0[2], P0[3], P0[4]*100, P0[5]*100, P0[6]*100)
  }
  x_spec <- spectra(x)[1,]
  if (is.null(transmittance_spectra))
  {
    t_spec <- NULL
  } else {
    t_spec <- spectra(transmittance_spectra)[1,]
  }
  res <- pracma::nelder_mead(x0 = P0, f = .inversionFUN, 
                             lb = c(10, 0, 0, -0.00000001, 0.005, 0.2),
                             ub = c(30, 100, 30, 100, 4, 1.8),
                             reflectance_spectra = x_spec, 
                             transmittance_spectra = t_spec, 
                             sam = sam, ...)
  res$xmin <- c(N = res$xmin[1]/10, Cab = res$xmin[2], Car = res$xmin[3], 
                Cbrown = res$xmin[4]/100, Cw = res$xmin[5]/100, Cm = res$xmin[6]/100)
  return(res)
}