##' Unit conversions
##'
##' Unit conversions
##' @title L2eV
##' @param wavelength wavelength in m
##' @return converted unit
##' @family conversion 
`L2eV` <- 
function(wavelength)
{
  with(constants, h * cel / ee / wavelength)
}

##' Unit conversions
##'
##' Unit conversions
##' @title eV2L
##' @param energy energy in eV
##' @family conversion 
`eV2L` <- function(energy)
{
  with(constants, h * cel / ee / energy)
}

##' Unit conversions
##'
##' Unit conversions
##' @title L2w
##' @param wavelength wavelength in m
##' @family conversion 
`L2w` <- 
function(wavelength)
{
  with(constants, 2*pi * cel / wavelength)
}
##' Unit conversions
##'
##' Unit conversions
##' @title t2eV
##' @param time time in s
##' @family conversion 
t2eV <- function(time){
  with(constants, h / (pi*time*ee))
}

