NULL
#' 
#' The water table recharge: the response unit  
#'
#' @param t time coordinate
#' @param d depth of unsaturated zone along the slope-normal direction
#' @param H soil depth
#' @param D soil water diffusivity
#' @param m maximum limit of summary truncation. Default is 100.
#' 
#' 
#' @note This function calcletes the water-table recharge rate in a hillslope assuming:
#'             
#' 1. Richards' Equation is linearized and reduced to the form of heat equation;
#' 
#' 2. The diffusion water-table rate is connectedwith soil pressure head according with eq. 13 (Cordano and Rigon, 2008);
#' 
#' @references
#' 
#' Cordano, E., and R. Rigon (2008), A perturbative view on the subsurface water pressure response at hillslope scale, Water Resour. Res., 44, W05407, doi:10.1029/2006WR005740.
#' \url{http://onlinelibrary.wiley.com/doi/10.1029/2006WR005740/pdf}
#' 
#' @export
#' 
#' @examples 
#' 
#' library(soilwater)
#' 
#' 
#' t <- seq(0,2,by=0.001)
#' d <- c(1,0.75,0.5,0.25)
#' val1 <- unitResponse(t, d = d[1], D = 1, H = 1, m = 500)
#' 
#' val2 <- unitResponse(t, d = d[2], D = 1, H = 1, m = 500)
#' 
#' val3 <- unitResponse(t, d = d[3], D = 1, H = 1, m = 500)
#' 
#' val4 <- unitResponse(t, d = d[4], D = 1, H = 1, m = 500)
#' 
#' 
#' 



unitResponse <- function(t,d=1,D=1,H=d,m=100) {
	
	sum <- 0 
	
	d <-d/H 
	
	tscale <- H^2/D 
	
	t <- t/tscale
	vect <- -m:m
	
	for (m in vect) {
		
		val <- (pi*t)^(-0.5)*(-1/(2*t)+(2*m*d+d)^2/(2*t)^2)*exp(-(2*m*d+d)^2/(4*t))
		sum <- sum+val
	}
	
	out <- sum*tscale 
	
	#
    # INPUT is DELTA rescaled with hydraulic conductivy
	# OUTPUT is rescaled with (d/tscale)
    #
	return(out)
	
	
}