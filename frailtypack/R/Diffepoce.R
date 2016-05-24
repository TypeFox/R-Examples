
Diffepoce <- function(epoce1, epoce2){

	if (class(epoce1)!="epoce") stop("Diffepoce allows only arguments of classe 'epoce'")
	if (class(epoce2)!="epoce") stop("Diffepoce allows only arguments of classe 'epoce'")
	if (!(all.equal(epoce1$pred.times,epoce2$pred.times))) stop("The two epoce objects should have the same prediction times")
	if (!(all.equal(epoce1$data,epoce2$data))) stop("The two epoce objects should be derived from the same dataset")

	DEPOCE <- as.double(epoce1$cvpol)-as.double(epoce2$cvpol)
	if (epoce1$new.data==TRUE){
		DEPOCE <- as.double(epoce1$mpol)-as.double(epoce2$mpol)
	}
	
	diff <- epoce1$IndivContrib-epoce2$IndivContrib
	diff_sq <- (epoce1$IndivContrib-epoce2$IndivContrib)^2

	DMPOL <- apply(diff,FUN=sum,MARGIN=2)/apply(epoce1$IndivContrib!=0,FUN=sum,MARGIN=2)

	omega <- (apply(diff_sq,FUN=sum,MARGIN=2)/apply(epoce1$IndivContrib!=0,FUN=sum,MARGIN=2)) - (DMPOL)^2
	omega <- omega/apply(epoce1$IndivContrib!=0,FUN=sum,MARGIN=2)

	TIinf <- as.vector(DEPOCE) - qnorm(0.975) * sqrt(omega)
	TIsup <- as.vector(DEPOCE) + qnorm(0.975) * sqrt(omega)

	out <- NULL
	out$new.data <- epoce1$new.data
	out$pred.times <- epoce1$pred.times
	out$DEPOCE <- DEPOCE
	out$TIinf <- TIinf
	out$TIsup <- TIsup

	class(out) <- c("Diffepoce")
	out
}