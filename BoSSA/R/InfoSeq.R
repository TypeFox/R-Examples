`InfoSeq` <-
function(X)
{
	paste(strsplit(X,"|",fixed=TRUE)[[1]][2],strsplit(X,"|",fixed=TRUE)[[1]][4],strsplit(X,"|",fixed=TRUE)[[1]][5],sep=",")
}

