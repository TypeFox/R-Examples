longX2 <- function(Freqs){
	# this is appropriate whenever P1 = 1 and P2 = 0 for all L markers
	if(class(Freqs)!="numeric") stop ("Freqs must be a vector of allele frequencies")
	if(sum(is.na(Freqs))>0) stop ("I can't handle NA's")
	if(sum(Freqs==0 | Freqs==1)>0) cat("Warning: frequencies of zero or one are present, causing division by zero. You should estimate Freqs as (x+1)/(n+2).\n")
	L <- length(Freqs) # number of markers
	M <- mean(Freqs,na.rm=TRUE) # MLE estimate of the proportionate contribution of P1
	V <- matrix(0,L,L)
	diag(V) <- M*(1-M)
	MSE <- t(Freqs-M)%*%solve(V)%*%(Freqs-M)/(L-1)
	VMi <- Freqs*(1-Freqs)*MSE
	X2 <- sum((Freqs-M)^2/VMi)
	list(test = data.frame(chisq=X2,df=L-1,p.val=pchisq(X2,df=L-1,lower.tail=FALSE)), chisq.res = (Freqs-M)/sqrt(VMi))
}

