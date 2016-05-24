#' Multifactor is a datastructure for survey items in the form of 'check all that apply'. Every response has multiple values.
#' @param values a list of vectors with response values
#' @param levels a vector with possible values 
#' @param labels a vector with labels
#' @param ordered ordered or unordered factor
#' @return a multifactor object
#' @aliases as.vector.multifactor expand.multifactor facdim is.multifactor as.vector.multifactor rep.multifactor [.multifactor [[.multifactor
#' @export
multifactor <- function(values, levels = unique(unlist(values)), labels=levels, ordered=TRUE){
	
	if(!all(na.omit(unlist(values)) %in% levels)){
		stop("Some values were not found in 'levels' at multifactor conversion.")
	}
	
	if(length(labels) != length(unique(labels))){
		warning("It seems like your multifactor as duplicate levels. That's not good:", paste(names(which(table(labels) > 1)), collapse=", "));
	}
	
	newvalues <- sapply(values, match, levels);
	newvalues <- sapply(newvalues, paste, collapse="+");
	newvalues[is.na(values)] <- NA;
	newlevels <- 1:length(levels);
	
	attr(newvalues, "levels") <- newlevels;
	attr(newvalues, "labels") <- labels;
	attr(newvalues, "ordered") <- ordered;
	class(newvalues) <- c("multifactor", "character");
	return(newvalues);	
}

#' @export
as.vector.multifactor <- function(x, mode){
	myvec <- unlist(strsplit(x, "+", fixed=TRUE));
	myfactor <- factor(myvec, attr(x, "levels"), attr(x, "labels"), attr(x, "ordered"));
	return(myfactor);
}

levelinfactor <- function(mylist, mylevel){
	if(length(mylevel) != 1) {
		stop("level has to be of length 1.");
	}
	return(sapply(sapply(mylist, "==", mylevel), any));
}

#' @export
expand.multifactor <- function(x){
	newvalues <- strsplit(x, "+", fixed=TRUE);
	mydf <- as.data.frame(sapply(attr(x, "levels"), levelinfactor, mylist=newvalues));
	colnames(mydf) <- attr(x, "labels");
	return(mydf);
}

#' @export
facdim <- function(x){
	return(sapply(strsplit(x, "+", fixed=TRUE), length));
}


#' @export
is.multifactor <- function(x){
	if("multifactor" %in% class(x)) {
		return(TRUE)
	} else {
		return(FALSE);
	}
}

#' @export
"[.multifactor" <- function(x, ..., drop = FALSE){
	y <- NextMethod("[")
	attr(y, "labels") <- attr(x, "labels")
	attr(y, "levels") <- attr(x, "levels")
	attr(y, "ordered") <- attr(x, "ordered")
	attr(y, "prompt_type") <- attr(x, "prompt_type")
	class(y) <- oldClass(x)
	lev <- levels(x)
	if (drop) 
		factor(y, exclude = if (any(is.na(levels(x)))) 
			NULL
		else NA)
	else y	
}

#' @export
"[[.multifactor" <- function (x, ...) {
	y <- NextMethod("[[")
	attr(y, "labels") <- attr(x, "labels")
	attr(y, "levels") <- attr(x, "levels")
	attr(y, "ordered") <- attr(x, "ordered")
	attr(y, "prompt_type") <- attr(x, "prompt_type")
	class(y) <- oldClass(x)
	y
}

#' @export
rep.multifactor <- function (x, ...) 
{
	y <- NextMethod()
	structure(y, class = class(x), levels = attr(x, "levels"), labels = attr(x, "labels"), ordered=attr(x, "ordered"), prompt_type = attr(x, "prompt_type"));
}


