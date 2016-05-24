syrjala0 <- function(coords, var1, var2, nsim, R=FALSE){
	
	datanames <- c(deparse(substitute(var1)),  deparse(substitute(var2)))
	
	if(R!=FALSE) {
		ppp1 = haz.ppp(cbind(coords,var1))
		ppp2 = haz.ppp(cbind(coords,var2))
		resultado <- syrjala.test(ppp1, ppp2, nsim)} else resultado <- syrjala(coords=coords, var1 =var1, var2=var2, nperm=nsim) 
		
		resultado$datanames <- datanames
		
	return(resultado)
}

