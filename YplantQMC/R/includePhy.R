
#'@rdname setPhy
#'@export includePhy
includePhy <- function(object,...){
	UseMethod("includePhy")
}

#'@rdname setPhy
#'@method includePhy plant3d
#'@S3method includePhy plant3d
includePhy.plant3d <- function(object, ...){
	phy <- setPhy(...)
	object$phy <- phy
return(object)
}


#'@rdname setPhy
#'@method includePhy plant3dlist
#'@S3method includePhy plant3dlist
includePhy.plant3dlist <- function(object, phydfr, ...){

	pfiles <- attributes(object)$pfiles
	nplants <- attributes(object)$nplants
	
	if(!is.data.frame(phydfr)){
		if(is.character(phydfr) && file.exists(phydfr))
			phydfr <- read.csv(phydfr)
		else
			stop("Provide 'phydfr' as a dataframe or name of CSV file")
	}	
	
	if(!"pfile" %in% names(phydfr))
		stop("phydfr must have a variable 'pfile' to match parameters to plants!")
	if(!"leafmodel" %in% names(phydfr))
		stop("phydfr must have a variable 'leafmodel' that specifies the leaf gas exchange model!")
	
	# Order of plants' pfiles in phydfr.
	porder <- match(pfiles, phydfr$pfile)
	if(any(is.na(porder)))
		stop("Some pfiles not in 'phydfr'. Please check inputs!")
	if(any(duplicated(porder)))
		stop("Some pfiles appear more than once in 'phydfr'. Please check inputs!")
	
	# names of parameters.
	parnames <- setdiff(names(phydfr), c("pfile","leafmodel"))
	
	# include.
	for(i in 1:nplants){
		p <- as.list(phydfr[porder[i],])
		p$pfile <- NULL
		l <- list(object=object[[i]], leafmodel=p$leafmodel, leafpars=p)
		l$leafpars$leafmodel <- NULL
		object[[i]] <- do.call("includePhy.plant3d", l)
	}

return(object)
}
