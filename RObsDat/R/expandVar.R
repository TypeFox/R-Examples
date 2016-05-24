expandVar <- function(var, nrow, ncol, checkID=FALSE, table=NULL){
	if(checkID){
		stopifnot(!is.null(table))
		theID  <- getID(table,var)
		dim(theID) <- dim(var)
	} else {
		theID <- var
	}
	if(is.matrix(theID)){
		nrow.var <- NROW(theID)
		ncol.var <- NCOL(theID)
		oneDim <- FALSE
	} else {
		nrow.var <- ncol.var <- length(theID)
		oneDim <- TRUE
	}
	#Same dimensions are required as were passed
	if(nrow.var == nrow & ncol.var == ncol & !oneDim){
		#make certain theID is an matrix
                dim(theID) <- c(nrow,ncol)
		return(theID)
	}

	#Different dimensions are required then  were passed

	if(ncol.var == 1 & nrow.var==1){
		return(matrix(theID, ncol=ncol, nrow=nrow))
	} else if(!is.matrix(theID)){
		if(nrow==ncol) warning(paste("Number of columns and rows equivalent for table'",table,"': unclear how to expand data! Assuming columns to match. To resolve this, pass arguments as matizes with the same dimension as your data.", sep=""))
		if(ncol==ncol.var){
			return(matrix(theID, ncol=ncol, nrow=nrow, byrow=TRUE))
		} else if(nrow==nrow.var){
			return(matrix(theID, ncol=ncol, nrow=nrow))
		} else {
			cat("Unexpected expansion of variable. Please implement\n")
			browser()
		}
	} else {
			cat("Unexpected expansion of variable. Please implement\n")
			browser()
	}
	stop("Should not be here")

}
