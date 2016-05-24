changeSupport <- function(from,to,min.coverage,FUN=NULL,weights.arg=NULL,
			  split.from=FALSE,merge.from=TRUE, ...)
UseMethod ('changeSupport')

setGeneric (name='changeSupport', def=function(from, to, min.coverage,
					       FUN=NULL, weights.arg=NULL,
					       split.from=FALSE, merge.from=TRUE,
					       ...)
	    standardGeneric('changeSupport') )
 
setMethod ('changeSupport',
	   signature(from='TimeIntervalDataFrame', to='TimeIntervalDataFrame',
		     min.coverage='numeric', FUN='ANY', weights.arg='ANY',
		     split.from='ANY', merge.from='ANY'),
	   definition=function (from, to, min.coverage, FUN=NULL, weights.arg=NULL,
				split.from=FALSE, merge.from=TRUE, ...)
{
	fun.args <- list (...)
	if (is.null (FUN) ) {
		if (length (fun.args) != 0 | !is.null(weights.arg) )
			warning ('Arguments passed to FUN are set to default.')
		if (homogeneous (from) ) {
			FUN <- mean
			fun.args <- list(na.rm = TRUE)
		} else {
			FUN <- weighted.mean
			weights.arg <- 'w'
			fun.args <- list(na.rm = TRUE)
		}
	}

	do.call(tapply, c(X=from, INDEX=to, FUN=FUN, fun.args, min.coverage=min.coverage,
		weights.arg=weights.arg, merge.X=merge.from, split.X=split.from,
		keep.INDEX=TRUE, simplify=TRUE))

} )

setMethod ('changeSupport',
	   signature(from='TimeIntervalDataFrame', to='character',
		     min.coverage='numeric', FUN='ANY', weights.arg='ANY', 
		     split.from='ANY', merge.from='ANY'),
	   definition= function (from, to, min.coverage, FUN=NULL,
				 weights.arg=NULL, split.from=FALSE,
				 merge.from=TRUE, ...)
{
	period <- POSIXctp(unit=to)
	return (changeSupport (from, period, min.coverage, FUN, weights.arg,
			       split.from, merge.from, ...) )
} )

setMethod ('changeSupport',
	   signature(from='TimeIntervalDataFrame', to='POSIXctp',
		     min.coverage='numeric', FUN='ANY', weights.arg='ANY', 
		     split.from='ANY', merge.from='ANY'),
	   definition= function (from, to, min.coverage, FUN=NULL,
				 weights.arg=NULL, split.from=FALSE,
				 merge.from=TRUE, ...)
{
	# if from and to have same base, no calculus to do 

	if (homogeneous (from) &&
	    continuous(from)  &&
	    to == period (from) )
		return (from)

	s <- min(start( from ))
	e <- max(end( from ))
	tzone <- timezone( from )
	to <- TimeIntervalDataFrame(s, e, period=to, timezone=tzone)
	
	return( changeSupport(from, to, min.coverage, FUN, weights.arg,
	      		      split.from, merge.from, ...) )
} )

