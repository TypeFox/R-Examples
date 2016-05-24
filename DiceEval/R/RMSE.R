RMSE <- function(Y,Ypred)
{
	Ypred <- as.numeric(Ypred)
	Y    <- as.numeric(Y)
	if(length(Ypred)!=length(Y)){
		stop("The entries must have the same length.")
	}
	return(sqrt(mean((Y-Ypred)^2)))
}