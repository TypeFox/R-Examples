#' Corrects a correlation for unreliability
#' 
#' \code{rel.correction} corrects obtained correlation for unreliability.
#' 
#' @param rxyi the obtained restricted correlation coefficient (obtained from
#' the incumbents).
#' @param rxxi the obtained reliability of x in the sample
#' @param ryyi the obtained reliability of y in the sample
#' @return returns a scalar that is the corrected correlation.
#' @author Dustin Fife
#' @export
#' @examples
#' 
#' 	rel.correction(rxyi=.3, rxxi=.8)
#' 	rel.correction(rxyi=.3, rxxi=.8, ryyi=.6)
#' 	rel.correction(rxyi=.3, ryyi=.6)		
rel.correction = function(rxyi, rxxi=NULL, ryyi=NULL){
	if (is.null(rxxi)) rxxi=1
	if (is.null(ryyi)) ryyi=1	
	rtp = rxyi/(sqrt(rxxi)*sqrt(ryyi))
	return(rtp)
}