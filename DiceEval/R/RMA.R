RMA <- function(Y,Ypred)
{
	Ypred <- as.numeric(Ypred)
	Y    <- as.numeric(Y)
	if(length(Ypred)!=length(Y)){
		stop("The entries must have the same length.")
	}
	tmp_ <- sort(abs(Y-Ypred)/sd(Y),index.return=TRUE,decreasing = TRUE)
	return(list(max.value=tmp_$x[1],max.data=tmp_$ix[1],index=tmp_$ix,error = tmp_$x))
}