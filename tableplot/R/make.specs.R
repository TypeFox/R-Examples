## make a collection of cell.specs for tableplot
#
#  each argument is repeated up to the required length
#  # of patterns = n (if specified), or the maximum length of any argument

#  list of lists to data.frame and back:
#    Given pats as a list of lists: 
#      pats.df <- do.call(rbind, lapply(pats, data.frame))
#      pats.df <- do.call(rbind.data.frame, pats)    ## works best
#      pats.df <- do.call(rbind, lapply(pats, data.frame))  ## works same -- real data frame
#      pats.df <- data.frame(do.call(rbind,lapply(pats,function(x) t(as.matrix(x,ncol=10)))))
#    Given pats.df: pats <- lapply(seq(along = rownames(pats.df)),function(i) as.list(pats.df[i, ]))

make.specs <- function(
  n=NULL,
  as.data.frame=FALSE,
  subset,
  ...
  ) {

	dots <- list(...)
	
	# Get the arguments and default values for cellgram()
	cell.args <- formals(cellgram)[-1]     # exclude cell value
	cell.arg.names <- names(cell.args)

	# If only a subset, select them
	# TODO:  what if selected names don't include all of dots??
	if (!missing(subset)) {
		wanted <- match(subset, names(cell.args))
#		wanted <- match(c(subset, names(dots)), names(cell.args))
		cell.args <- cell.args[!is.na(wanted)]
	}
	cell.arg.names <- names(cell.args)
	nargs <- length(cell.args)

	# find max length among dots arguments
	if(!is.null(n)) len <- n
	else {
		len <- max(unlist( lapply( dots, FUN=length )))
	}
	
	# check for dots arguments not among arguments to cellgram
	if (! all(names(dots) %in% names(cell.args)))
		warning(paste(names(dots)[!names(dots) %in% names(cell.args)], "are not cellgram arguments; ignored."))
		
#  # replicate the elements in each, to required length
	args.in.dots <- match(names(cell.args), names(dots))
	for (i in seq_along(cell.args)) {
		if (!is.na(args.in.dots[i]))
			cell.args[[i]] <- rep(dots[[args.in.dots[i]]], length.out=len)
		else cell.args[[i]] <- rep(cell.args[[i]], length.out=len)
		}

#browser()

  result <- as.data.frame(cell.args, stringsAsFactors=FALSE)
  if (as.data.frame) return(result)
	# reshape as list of lists
  result <- lapply(split(result, 1:nrow(result)), as.list)
  return(result)

}