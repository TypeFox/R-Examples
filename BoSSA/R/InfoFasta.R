`InfoFasta` <-
function(X,Y)
{
	nom <- names(X)
	info <- lapply(nom,InfoSeq)
	info <- as.matrix(info)
	write.table(info,paste(substring(as.character(Y),1,5),".txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
}

