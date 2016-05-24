

roc <- function(x, ...) UseMethod("roc")



roc.integer <- function(x, ...) {
  roc.numeric(x, ...)
}



roc.numeric <- function(x,labels,labpos,thres=NULL,...){
	id <- labels == labpos
	if(sum(id) < 1)
		stop("\n wrong level specified !\n")
	if(is.null(thres))
		thres <- sort(unique(x))
	else
		thres <- sort(unique(thres))
	ans <- .Call("rocDaim",as.numeric(x),
				as.numeric(id),
				as.numeric(thres),
				as.numeric(0),
				PACKAGE="Daim")
	class(ans) <- "Daim.vector"
	ans
}



roc.matrix <- function(x,labels,labpos,thres=NULL,...){
	id <- labels == labpos
	if(sum(id) < 1)
		stop("\n wrong level specified !\n")
	if(!is.null(thres))
		thres <- sort(unique(thres))
	ans <- apply(x,2,
				function(y){
				if(is.null(thres))
					thres <- sort(unique(y))
				erg <- .Call("rocDaim",as.numeric(y),
						as.numeric(id),
						as.numeric(thres),
						as.numeric(0),
						PACKAGE="Daim")
				}
			)
	class(ans) <- "Daim.list"
	ans
}



roc.data.frame <- function(x, ...) {
  roc.matrix(x, ...)
}



roc.default <- function(x, ...) {
  stop(paste("Do not know how to handle objects of class", class(x)))
}








