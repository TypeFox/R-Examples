"timedep.names" <- function(x)
{
	pos1 <- grep("p",unlist(strsplit(x,split="")))[1]+2
	pos2 <- grep(")",unlist(strsplit(x,split="")))[1]-1
	pos3 <- grep(")",unlist(strsplit(x,split="")))[1]+1
	pos4 <- length(unlist(strsplit(x,split="")))
	return(paste(substr(x,start=pos1,stop=pos2),substr(x,start=pos3,stop=pos4),sep=""))
}