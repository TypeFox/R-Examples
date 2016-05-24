bootvar <-
function(X.quanti=NULL,X.quali=NULL)
{
	if (!is.null(X.quanti)&& !is.null(X.quali))
	{
		n <- nrow(X.quanti)
		indice <- sample(1:n,n,replace=TRUE)
		Xboot.quanti <- X.quanti[indice,]
		Xboot.quali <- X.quali[indice,]
	}
	if (!is.null(X.quanti)&& is.null(X.quali))
	{
		n <- nrow(X.quanti)
		indice <- sample(1:n,n,replace=TRUE)
		Xboot.quanti <- X.quanti[indice,]
		Xboot.quali <- NULL
	}
	if (is.null(X.quanti)&& !is.null(X.quali))
	{
		n <- nrow(X.quali)
		indice <- sample(1:n,n,replace=TRUE)
		Xboot.quali <- X.quali[indice,]
		Xboot.quanti <- NULL
	}

	list(Xboot.quanti=Xboot.quanti,Xboot.quali=Xboot.quali)
}

