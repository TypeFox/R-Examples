#' Merge two ffdf by common columns, or do other versions of database join operations. 
#'
#' Merge two ffdf by common columns, or do other versions of database join operations.
#' This method is similar to \code{\link[base]{merge}} in the \code{base} package but only allows inner and left outer joins.
#' Note that joining is done based on \code{\link{ffmatch}} or \code{\link{ffdfmatch}}: only the first element
#' in \code{y} will be added to \code{x}; and since \code{\link{ffdfmatch}} works by \code{\link[base]{paste}}-ing together a key, 
#' this might not be suited if your key contains columns of vmode double.
#'
#' If a left outer join is performed and no matching record in x is found in y, columns with vmodes 'boolean', 'quad', 'nibble', 'ubyte', 'ushort' are
#' coerced respectively to vmode 'logical', 'byte', 'byte', 'short', 'integer' to allow NA values.
#'
#' @method merge ffdf
#' @export
#' @export merge.ffdf
#' @example ../examples/merge.R
#' @param x an ffdf
#' @param y an ffdf
#' @param by specifications of the common columns. Columns can be specified by name, number or by a logical vector.
#' @param by.x specifications of the common columns of the x ffdf, overruling the by parameter
#' @param by.y specifications of the common columns of the y ffdf, overruling the by parameter
#' @param all see \code{\link{merge}} in R base
#' @param all.x if TRUE, then extra rows will be added to the output, one for each row in x that has no matching row in y. 
#' These rows will have NAs in those columns that are usually filled with values from y. The default is FALSE, so that only rows with data from both x and y are included in the output. 
#' @param all.y similar as all.x
#' @param sort logical, currently not used yet, defaults to FALSE.
#' @param suffixes character(2) specifying the suffixes to be used for making non-by names() unique.
#' @param incomparables values which cannot be matched. See \code{\link{match}}. Currently not used.
#' @param trace logical indicating to show on which chunk the function is computing
#' @param ... other options passed on to \code{\link[ff]{ffdfindexget}}
#' @return 
#' an ffdf
#' @seealso \code{\link{merge}}
merge.ffdf <- function(x, y, by = intersect(names(x), names(y)), by.x = by, by.y = by, all=FALSE, all.x = all, all.y = all, sort=FALSE, 
	suffixes = c(".x",".y"), incomparables=NULL, trace=FALSE, ...){
	if (!is.null(incomparables)){
    .NotYetUsed("incomparables != NULL")
  }
  if((all.x == TRUE & all.y == TRUE) | (all.y == TRUE & all.x == TRUE)){
  	stop("merge.ffdf only allows inner joins")
  }
  if (!all(by.x %in% names(x))){
    stop("'",by.x[!by.x %in% names(x)], "' is missing from x")
  }
	if (!all(by.y %in% names(y))){
	  stop("'", by.y[!by.y %in% names(y)], "' is missing from y")
	}
	## First find the first match where the left hand side ffdf is in the right hand side ffdf
  if(length(by.x) == 1){
  	matchidx <- eval(ffmatch(x[[by.x]], y[[by.y]], trace = trace))
  }else{
  	matchidx <- eval(ffdfmatch(x[by.x], y[by.y], trace = trace))
  }  
  updatepositions <- ffwhich(matchidx, !is.na(matchidx))
  ## Define the type of join
  if(all.x == FALSE & all.y == FALSE){  	
	  jointype <- "inner"
	  if(length(updatepositions) == 0){
	  	## No match is found and an inner join is requested
	  	warning("No match found, returning NULL as ffdf can not contain 0 rows")
	  	return(NULL)
	  }
  }else if(all.x == TRUE & all.y == FALSE){
    jointype <- "left"
  	if(any(is.na(matchidx)) == FALSE){
  		jointype <- "inner"
  	}
  }
  xvars <- colnames(x)
  newxvars <- ifelse(!xvars %in% by.x & xvars %in% colnames(y), paste(xvars, suffixes[1], sep=""), xvars)
  yvars <- colnames(y)[!colnames(y) %in% by.y]
  newyvars <- ifelse(yvars %in% colnames(x), paste(yvars, suffixes[2], sep=""), yvars)
  if(length(yvars) == 0){
  	stop("merge.ffdf requires at least 1 column in y not in by.y")
  }
  ## JOIN  
  if(jointype == "inner"){
  	##
  	## Everything from the left is matched, we don't need to change the vmode  	
  	if(trace) {
    	message(sprintf("%s, found match indexes, now starting to add y to x", Sys.time()))
  	}	
  	res <- ffbaseffdfindexget(y[yvars], index=matchidx[updatepositions], ...)
  	## cbind it to the existing x ffdf if there is a match
  	x <- ffbaseffdfindexget(x, index=updatepositions, ...)
  	colnames(x) <- newxvars
  	for(i in 1:length(yvars)){
  		x[[newyvars[i]]] <- res[[yvars[i]]]
  	}	  	 
  }else if(jointype == "left"){
  	##
  	## We need to fix the vmode if needed as we are adding NA's
  	if(trace){
    	message(sprintf("%s, found match indexes, now starting to add y to x and coercing if needed", Sys.time()))
  	}  	
  	tocoerce <- coerce_to_allowNA(vmode(y[yvars]))  	
  	## First set everything to NA  	
  	for(i in 1:length(yvars)){
  		x[[newyvars[i]]] <- ff::clone(y[[yvars[i]]], initdata=NA, length = nrow(x), vmode = tocoerce$coerceto[[yvars[i]]])
  	} 
  	## Next coerce the right hand side if needed to allow NA values
  	tocoerce.yvars <- yvars[tocoerce$x != tocoerce$coerceto]
  	if(length(tocoerce.yvars) > 0){
  		warning(paste("coercing column ", paste(tocoerce.yvars, collapse=", "), " to a higher vmode to allow NA's"), sep="")
  		for(measure in tocoerce.yvars){
 				y[[measure]] <- ff::clone(y[[measure]], vmode = tocoerce$coerceto[[measure]])
  		}
  	} 		
  	## Next update the existing x ffdf with the joined value if there is a match  	
  	if(length(updatepositions) > 0){
	  	matchidx.nonmissing <- matchidx[updatepositions]
  		x[updatepositions, newyvars] <- ffbaseffdfindexget(y[yvars], index=matchidx.nonmissing)
  	}  	
  }    
	x
}

