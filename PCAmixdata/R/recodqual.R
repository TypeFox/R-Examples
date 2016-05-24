recodqual <-
function(X,rename.level=FALSE)
	{
		X <- as.matrix(X)
		GNA <- tab.disjonctif.NA(X,rename.level)
		G <- replace(GNA,is.na(GNA),0)
		ns <- apply(G,2,sum)
		nmiss <- apply((apply(GNA,2,is.na)),2,sum)
		n <- nrow(X)
		if(sum((n-nmiss)==ns)!=0) stop("There are columns in X.quali where all the categories are identical")
		return(G)	
	}

