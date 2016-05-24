R2 <- function(Y,Ypred)
{
	Ypred <- as.numeric(Ypred)
	Y     <- as.numeric(Y)
	if(length(Ypred)!=length(Y)){
		stop("The entries must have the same length.")
	}
	return(1 - mean((Y-Ypred)^2)/mean((Y-mean(Y))^2))
}