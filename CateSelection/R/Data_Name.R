Data_Name <-
function(dat,heter=NULL,...){
	if(is.null(heter))heter <- FALSE
	ColNum <- length(dat)
	datColNames <- colnames(dat)   ##save original column names
	datColNames <- datColNames[-ColNum]
	datColNames <- t(t(datColNames))
	ColNames <- numeric() 
	for(i in 1:(ColNum-1)){ColNames[i] <- paste("x",i,sep="")}  ##change column names
	ColNames[ColNum] <- "y"
	colnames(dat) <- ColNames
	if(heter)for(i in 1:(length(dat)-1))dat[,i] <- replace(dat[,i],which(dat[,i]=="H"),NA) ##relace "H" to NA
	results <- list(datColNames=datColNames,dat=dat)
	return(results)
}
