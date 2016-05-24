`rename` <-
function(X)
{
	nom <- names(X)
	nom2 <- lapply(nom,AccNumb)
	names(X) <- nom2
	X
}

