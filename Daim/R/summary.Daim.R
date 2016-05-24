
summary.Daim <- function(object,...)
{
	meth <- unlist(object$method)
	meth.names <- names(meth)
	meth <- as.vector(meth)
	auc <- auc(object)
	if(class(object)[2] != "cv"){
		ans <- list(call=object$formula,
					method=list(meth.names=meth.names,meth=meth,cutoff=object$cutoff),
					error=data.frame(err.632p=object$err632p, err.632=object$err632,
									 err.loob=object$errloob, err.app=object$errapp),
					sensitivity=data.frame(sens.632p=object$sens632p, sens.632=object$sens632,
										   sens.loob=object$sensloob, sens.app=object$sensapp),
					specificity=data.frame(spec.632p=object$spec632p, spec.632=object$spec632,
										   spec.loob=object$specloob, spec.app=object$specapp),
					AUC=data.frame(auc[1:4]))
	}
	else{
		ans <- list(call = object$formula,
					method = list(meth.names=meth.names, meth=meth, cutoff=object$cutoff),
					error = data.frame(err.loob=object$errloob, err.app=object$errapp),
					sensitivity = data.frame(sens.loob=object$sensloob, sens.app=object$sensapp),
					specificity = data.frame(spec.loob=object$specloob, spec.app=object$specapp),
					AUC = data.frame(auc[1:2]))
	}
	class(ans) <- c("summary.Daim",class(object)[2])
	ans
}



print.summary.Daim <- function(x, digits=max(3, getOption("digits") - 3), ...)
{
	method <- paste(paste(x$method$meth.names,x$method$meth,sep=" = "),"",sep=",")
	method[length(method)] <- gsub(",",".",method[length(method)])
	cat("\nPerformance of the classification obtained by:\n")
	cat("\nCall:\n ")
	print(formula(x$call))
	cat("\nDaim parameters: \n")
	cat(method,labels=" ",fill=TRUE)
	cat("\nResult: \n")
	error <- format(round(x$error,digits=digits),nsmall=digits)
	sens <- format(round(x$sensitivity,digits=digits),nsmall=digits)
	spec <- format(round(x$specificity,digits=digits),nsmall=digits)
	AUC <- format(round(x$AUC,digits=digits),nsmall=digits)
	if(class(x)[2] != "cv"){
		cat("---------------------------------------------------------------\n")
		cat("| Method:       ","|",paste(".632+ ",".632  ","loob  ","apparent ",sep="  |  "),"|\n")
		cat("===============================================================\n")
		cat("| Error:        ","|",paste(error[1],error[2],error[3],error[4],sep="  |  "),"   |\n")
		cat("---------------------------------------------------------------\n")
		cat("| Sensitivity:  ","|",paste(sens[1],sens[2],sens[3],sens[4],sep="  |  "),"   |\n")
		cat("---------------------------------------------------------------\n")
		cat("| Specificity:  ","|",paste(spec[1],spec[2],spec[3],spec[4],sep="  |  "),"   |\n")
		cat("---------------------------------------------------------------\n")
		cat("| AUC           ","|",paste(AUC[1],AUC[2],AUC[3],AUC[4],sep="  |  "),"   |\n")
		cat("---------------------------------------------------------------\n\n")
		invisible(x)
	}
	else{
		cat("-----------------------------------------\n")
		cat("| Method:       ","|",paste("  cv  ","apparent ",sep="  |  "),"|\n")
		cat("=========================================\n")
		cat("| Error:        ","|",paste(error[1],error[2],sep="  |  "),"   |\n")
		cat("-----------------------------------------\n")
		cat("| Sensitivity:  ","|",paste(sens[1],sens[2],sep="  |  "),"   |\n")
		cat("-----------------------------------------\n")
		cat("| Specificity:  ","|",paste(spec[1],spec[2],sep="  |  "),"   |\n")
		cat("-----------------------------------------\n")
		cat("| AUC           ","|",paste(AUC[1],AUC[2],sep="  |  "),"   |\n")
		cat("-----------------------------------------\n\n")
		invisible(x)
	}
}




summary.Daim.vector <- function(object, ...){
	best.id <- which.max((1-object$FPR)+object$TPR-1)
	best.cut <- object$cutoff[best.id]
	best.FPR <- object$FPR[best.id]
	best.TPR <- object$TPR[best.id]
	AUC <- auc(1-object$FPR,object$TPR)
	ans <- list(data=object, N=length(object$FPR), 
				best.cut=best.cut, FPR=best.FPR,
				TPR=best.TPR, AUC=AUC)
	class(ans) <- "summary.Daim.vector"
	ans
}




summary.Daim.list <- function(object, ...){
	ans <- lapply(object, summary.Daim.vector)
	class(ans) <- "summary.Daim.list"
	ans
}



print.summary.Daim.list <- function(x, digits=max(3, getOption("digits") - 3), ...)
{
	names.erg <- names(x)
	if(is.null(names.erg))
	names.erg <- paste("Col.",1:length(x),sep="")
	erg <- matrix(unlist(lapply(x,Daim.print.list)),nrow=5)
	dimnames(erg) <- list(rep("",5),names.erg)
	class(erg) <- "table"
	print(erg)
}



print.summary.Daim.vector <- function(x, digits=max(3, getOption("digits") - 3), ...)
{
	cat("\n Number of cut-points:",x$N-1,"\n")
	cat(" Best cut-point      :",format(round(x$best.cut,digits=digits),nsmall=digits),"\n")
	cat(" Sensitivity         :",format(round(x$TPR,digits=digits),nsmall=digits),"\n")
	cat(" Specificity         :",format(round(1-x$FPR,digits=digits),nsmall=digits),"\n")
	cat(" AUC                 :",format(round(x$AUC,digits=digits),nsmall=digits),"\n")
}



Daim.print.list <- function(x, digits=max(3, getOption("digits") - 3), ...)
{
	su.names <- c("Number of cut-points:","Best cut-point      :",
				  "Sensitivity         :", "Specificity         :",
				  "AUC                 :")
	erg <- paste(su.names,format(round(unlist(x[2:6]),digits=digits),nsmall=digits),"  ")
	class(erg) <- "table"
	erg
}



