lsdata <- function(fnam = '.Rdata') {
# cheap way to view objects in an Rdata file
#make CRANripley happy...
	x <- load(fnam, envir = environment())
	return(x)
#nobody likes 'attach' any more
	#datfoo<-attach(fnam)
	#foo<- ls(envir=datfoo)
	#detach(paste('file:',fnam,sep='',collapse=''),character.only=TRUE)
	#return(foo)
	}
