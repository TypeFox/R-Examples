inspect <-
function(x, what, retain = FALSE){
	z <- x$tables
	if(is.list(what)){
		those <- sapply(names(what), function(x) grep(x, names(z)))
		check <- sapply(those, length)
		if(length(which(check!=1))!=0){stop(paste("The items ", which(check!=1), "in 'what' are ambiguously defined. Please check your function call"))}
		for(i in c(1:length(what))){
			z[[those[i]]] <- z[[those[i]]][-what[[i]],]
		}
		if(retain){x$tables.orig <- x$tables}
		x$tables <- z
		return(x)
	}
	else{
		if(is.character(what)){what <- sapply(what, function(x) grep(x, names(z)))}
		if(length(what)==0){stop(paste("'what' is ambiguously defined. Please check your function call"))}
		return(z[what])
	}
}