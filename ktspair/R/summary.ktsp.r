summary.ktsp <- function(object,select=NULL,printall=FALSE,...){
	ktspobj <- object
	grp <- ktspobj$grp
	k <- ktspobj$k
	grplabels <- character(length(grp))
	grplabels[grp==0] <- ktspobj$labels[1]
	grplabels[grp==1] <- ktspobj$labels[2]

	if(!is.null(select) && (select <1 || select > k)){stop("The selected pair of genes is not available.")}

	if(!is.null(select)){
		cat(paste("Data for TSP:", select,"\n"))
		prediction <- predict(ktspobj, selec = select)
		print(table(prediction,grplabels,dnn=list(paste("1(Gene",rownames(ktspobj$ktspdat)[select]," < Gene", rownames(ktspobj$ktspdat)[(select + k)],")"),"Group Labels")))
	}

	if(is.null(select) & printall==TRUE){
		cat(paste("There are",k,"TSPs\n\n"))
		for(i in 1:k){
			cat(paste("Data for TSP:", i,"\n"))
			prediction <- predict(ktspobj, select = i)
			print(table(prediction,grplabels,dnn=list(paste("1(Gene",rownames(ktspobj$ktspdat)[i]," < Gene", rownames(ktspobj$ktspdat)[(i + k)],")"),"Group Labels")))
			cat("\n\n")
			readline("Hit return for next TSP")
		}
	}
	if(is.null(select) & printall==FALSE){
		label <- numeric(length(grp))
		cat(paste("There are ",k,"TSPs\n\n"))    
		prediction <- predict(ktspobj)
		cat(paste("Data for the k-TSP\n"))
		print(table(prediction,grplabels, dnn=list("Prediction Labels","Group Labels")))  
	}
}


