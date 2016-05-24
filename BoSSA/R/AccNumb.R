`AccNumb` <-
function(X)
{
	strsplit(strsplit(X,"|",fixed=TRUE)[[1]][4],".",fixed=TRUE)[[1]][1]
}

