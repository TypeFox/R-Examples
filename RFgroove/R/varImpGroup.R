varImpGroup <-
function(object, xdata, ngroups=length(nvarGroup), nvarGroup, idxGroup, groupsNames=names(nvarGroup), normalize=(length(unique(nvarGroup))!=1) ) {

	##
	## Version 5 - 07-07-14
	##

	# cat("normalize =", normalize, "\n")

	# For temporal importance
	if(is.list(idxGroup))
		idxGroup <- unlist(idxGroup)

	if(!is.numeric(idxGroup) | min(idxGroup)!=0)
		stop("'idxGroup' must contain the indexes of the grouped variables starting from 0")



	if(is.null(forest <- object$forest))
		stop("Error: keepForest")

	if(is.null(object$inbag))
		stop("Error: keepInbag")

	if(class(xdata)=="data.frame")
		xdata <- as.matrix(xdata)

	if(!is.null(groupsNames)){
		if(length(nvarGroup) != length(groupsNames))
			stop("Error: length nvarGroup and groupsNames")
	}
	
	M <- object$ntree
	Y <- object$y
	n <- nrow(xdata)
	p <- ncol(xdata)


	if(length(nvarGroup)!=ngroups){ #length(nvarGroup)>p | 
		# cat("length(nvarGroup) =", length(nvarGroup), " p =", p, " ngroups =", ngroups, "\n")
		stop("Error: number of groups")
	}


	xdata <- t(xdata)
	storage.mode(xdata) <- "double"

	if(is.factor(Y)){
		# cat("VarImpGroup for classification\n")
		obj <- .C(	"R_varImpGroup", varImpGroup=numeric(ngroups), xdata, as.integer(Y), n=as.integer(n), p=as.integer(p), ntree=as.integer(M), 
					as.integer(aperm(forest$treemap, c(2,1,3))), as.integer(forest$nodestatus), as.double(forest$xbestsplit),
					as.integer(forest$bestvar), as.integer(forest$nodepred), as.integer(forest$ndbigtree), as.integer(forest$ncat),
					as.integer(forest$maxcat), as.integer(length(unique(Y))), ngroups=as.integer(ngroups), as.integer(nvarGroup), as.integer(max(nvarGroup)),
					as.integer(idxGroup), as.integer(object$inbag), as.integer(forest$nrnodes), PACKAGE="RFgroove")
	}else{
		# cat("VarImpGroup for regression\n")
		obj <- .C(	"R_varImpGroup_Reg", varImpGroup=numeric(ngroups), xdata, as.double(Y), n=as.integer(n), p=as.integer(p), ntree=as.integer(M), 
					as.integer(forest$leftDaughter), as.integer(forest$rightDaughter), as.integer(forest$nodestatus), as.double(forest$xbestsplit),
					as.integer(forest$bestvar), as.double(forest$nodepred), as.integer(forest$ndbigtree), as.integer(forest$ncat), 
					as.integer(max(forest$ncat)), ngroups=as.integer(ngroups), as.integer(nvarGroup), as.integer(max(nvarGroup)), 
					as.integer(idxGroup), as.integer(object$inbag), as.integer(forest$nrnodes), PACKAGE="RFgroove")
	}

	if(normalize){
		IMP <- obj$varImpGroup / nvarGroup
	}else{
		IMP <- obj$varImpGroup
	}
	if(is.null(groupsNames)){
		names(IMP) <- paste("Gr", 1:ngroups, sep="")
	}else{
		names(IMP) <- groupsNames
	}
	class(IMP) <- "importance"
	return(IMP)
}
