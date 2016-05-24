###################################################################################
#' Calculates the first derivative of signal data
#'
#' @param DATA	 			The matrix of signals for which the derivative is calculated (one column per signal)
#' @param axType			Type of interpolation: "ORD1" (default) first order, "SCD" second ,"TRD" third, "ORD3a" cubic polynom
#'
#' @note setting axType to invalid values will lead to first order interpolation.
#'
#' @return matrix \code{DATA} \cr
#' - \code{DATA} contains the derivative signals, with the same structure as the input data.
#' @export
#' @importFrom stats filter
###################################################################################
sfaTimediff <- function (DATA, axType="ORD1"){
	#if (is.null(axType)){
	#	# first order interpolation
	#	DATA=filter(DATA,filter=c(1, -1))#matlab version DATA=filter([1 -1], [1], DATA);
	#}	
	#else{ 
		if (axType=="ORD3a"){
			# interpolation by cubic polynom
			DATA=filter(DATA,filter=c(2,-9, 18, -11))#DATA=filter([2 -9 18 -11], [1], DATA)./6;
		}
		else if (axType=="SCD"){
			DATA=filter(DATA,filter=c(1, -2, 1))#DATA=filter([1 -2 1], [1], DATA);
		}
		else if (axType=="TRD"){
			DATA=filter(DATA,filter=c(-1, 3, -3, 1))#DATA=filter([-1 3 -3 1], [1], DATA);
		}
		else{
			# first order interpolation
			DATA=filter(DATA,filter=c(1, -1))#DATA=filter([1 -1], [1], DATA);
		}
	#}
	#remove NA values (number of removed depends on filter order)
	if  (customSize(DATA)[1]>1){
		DATA=DATA[!is.na(DATA)[,1],]
	}else{
		DATA=DATA[!is.na(DATA)]
	}
	
  	return(DATA)
}