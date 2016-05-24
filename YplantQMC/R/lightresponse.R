#'Non-rectangular hyperbola
#'
#'@description A simple light response function that predicts leaf photosynthesis from
#'absorbed PAR.
#'
#'
#'@param PAR Photosynthetically active radiation (mu mol m-2 s-1)
#'@param Amax Maximum assimilation rate (the asymptote) (mu mol m-2 s-1)
#'@param phi Quantum yield (slope at PAR = 0) (mol mol-1)
#'@param theta Shape of light response curve (0 = rectangular hyperbola, 1 =
#''Blackman' response).
#'@param Rd Dark respiration (**positive** value).
#'@param \dots Further arguments are ignored.
#'@return Returns a dataframe with : \describe{ \item{A}{Net assimilation rate
#'(mu mol m-2 s-1)} }
#'@author Remko Duursma
#'@seealso \code{\link{setPhy}}
#'@keywords misc
#'@export
lightresponse <- function(PAR, Amax, phi, theta, Rd, ...){
                     

   A <- (phi*PAR+Amax-((phi*PAR+Amax)^2-4*theta*phi*PAR*Amax)^0.5)/(2*theta) - Rd
	   
return(data.frame(A=A))	   
}
