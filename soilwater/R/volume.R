# TODO: Add comment
# 
# Author: ecor
###############################################################################

NULL
#'
#' Water volume in function of water-table depth or height 'swc', Hydraulic Conductivity 'khy' , Soil Water Capacity 'cap' , Soil Water (Hydraulic) Diffusivity 'diffusivity'
#' 
#' 
#' @param d water-table depth (under surface) 
#' @param h water-table heigth (over bedrock)
#' @param H soil thickness
#' @param Gamma liner coefficient for hydrostatic profile (Default is 1)
#' @param nstep number of vertical spatial cells. Default is 100
#' @param soilwaterretentioncurve function describing the soil water retention curve. Default is  \code{\link{swc}}
#'
#' @param ... parametes for \code{soil.water.retention.curve}
#' 
#' @seealso \code{\link{swc}}
#' 
#' @note The water volume per topographical area unit obtained by vertical integration off soil water content profile 
#' 
#' 
#' @export




watervolume <- function(d=H-h,H=1,h=NA,nstep=100,Gamma=1,soilwaterretentioncurve=swc,...)  {
	
	profile <- seq(from=0,to=H,length.out=nstep)
	
	dz <- profile[-1]-profile[-nstep]
	z  <- (profile[-1]+profile[-nstep])/2
	
	psi <- (z-d)*Gamma
	
	theta <- soilwaterretentioncurve(psi,z=z,...)
	
	out <- t(theta) %*% dz
	
#	out <- NULL
	
	
	return(out)
	
}