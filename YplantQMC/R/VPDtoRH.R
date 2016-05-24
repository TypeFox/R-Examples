VPDtoRH <- function(VPD, TdegC, Pa=101){
	esatval <- esat(TdegC)
	e <- pmax(0, esatval - VPD*1000)
	RH <- 100 * e/esatval
return(RH)
}
