
print.Daim <- function(x, digits=max(3, getOption("digits") - 3), ...){
	meth <- unlist(x$method)
	names.meth <- names(meth)
	meth <- as.vector(meth)
	method <- paste(paste(names.meth,meth,sep=" = "),"",sep=",")
	method[length(method)] <- gsub(",",".",method[length(method)])
	cat("\nPerformance of the classification obtained by:\n")
	cat("\nCall:\n ")
	print(formula(x$formula))
	cat("\nDaim parameters: \n")
	cat(method,labels=" ",fill=TRUE)
	cat("\nResult: \n")
	if(class(x)[2] != "cv"){
		erg <- format(round(c(x$err632p,x$err632,x$errloob,x$errapp),digits=digits),nsmall=digits)
		cat("---------------------------------------------------------\n")
		cat(" Error:   ","|",paste(".632+ ",".632  ","loob  ","apparent ",sep="  |  "),"|\n")
		cat("           ----------------------------------------------\n")
		cat("          ","|",paste(erg[1],erg[2],erg[3],erg[4],sep="  |  "),"   |\n")
		cat("---------------------------------------------------------\n")
		invisible(x)
	}
	else{
		erg <- format(round(c(x$errloob,x$errapp),digits=digits),nsmall=digits)
		cat("-----------------------------------\n")
		cat(" Error:   ","|",paste(" cv   ","apparent ",sep="  |  "),"|\n")
		cat("           ------------------------\n")
		cat("          ","|",paste(erg[1],erg[2],sep="  |  "),"   |\n")
		cat("-----------------------------------\n")
		invisible(x)
	}
}



