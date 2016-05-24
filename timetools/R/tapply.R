# definition of the generic (S4) function
setGeneric (name='tapply')

# tapply method to apply FUN over TimeIntervalDataFrame splitted by 
# another TimeIntervalDataFrame
# this method is the main one : other tapply methods (over Time*DataFrame)
# call that one after an internal change of input datas.

setMethod (
	'tapply',
	signature(X='TimeIntervalDataFrame', INDEX='TimeIntervalDataFrame'),
	definition=function (X, INDEX, FUN, ..., min.coverage=1,
			     weights.arg=NULL, merge.X=TRUE,
			     split.X=FALSE, keep.INDEX=TRUE,
			     simplify=TRUE)
{
	# get the function and specified args to apply

	fun.args <- list( ... )
	FUN <- match.fun( FUN )

	# get useful time properties of the INDEX object

	s.INDEX <- start(INDEX)
	e.INDEX <- end(INDEX)

	# min.coverage must be between 0 and 1 or NA, validity test.

	if( (min.coverage < 0 | min.coverage > 1) & !is.na(min.coverage) )
		stop ("'min.coverage' must be between [0-1] or NA." )

	# X is split over INDEX using the split method defined
	# in the package.

	splitted <- split(X, INDEX, split.x=split.X, keep.f = FALSE)

	# if min.coverage is NA, no problem.
	# if it is valid and between 0 and 1, splitted data mustn't 
	# overlap : if this is the case coverage of data can not be 
	# calculated.

	if( !is.na( min.coverage ))
		if( any(sapply(splitted, overlapping)) )
			stop ("Overlapping data not allowed. Set min.coverage to NA.")

	# if merge.X is FALSE, splitted datas are not allowed to have
	# more than one line. When a splitted data has more than one line
	# all its values are set to NA. When the function is applied 
	# (next step), the aggregation result should be NA (sure ?)

	if( !merge.X )
		splitted[sapply(splitted, nrow) > 1] <-
			lapply(splitted[sapply(splitted, nrow) > 1],
			       function(x){x[TRUE] <- NA ; x})

	# the 'calc' subfunction is in charge of applying FUN and its args
	# to each splitted datas.
	# TODO : should this function be defined externally to the tapply FUN?

	calc <- function(x, s, e, mc, ws.a, FUN, f.as)
	{
		# x	: a TimeIntervalDataFrame (from the splitted datas)
		# s	: start of the resulting interval (get from INDEX)
		# e	: end of the resulting interval (get from INDEX)
		# mc	: min.coverage. Used to check that time coverage
		#	of 'x' is enough against time interval defined between 
		#	's' and 'e'.
		# ws.a	: if a 'weights.arg' is defined, its name, otherwise
		#	NULL
		# FUN	: the function to apply
		# f.as	: arguments of FUN sepcified by the user

		# if 'x' is empty, a TimeIntervalDataFrame beginning at 's'
		# ending at 'e' with all 'x' variables set to NA is returned.
		if( nrow(x) == 0 )
			return( merge(x, TimeIntervalDataFrame(s, e)) )

		# if necessary, weights for each row of 'x' are calculated
		if(!is.null( ws.a ) | !is.na( mc ))
			w <- as.numeric(end(x)) - as.numeric(start(x))
		if(!is.null( ws.a ))
			f.as[[ws.a]] <- w

		# duration of the returned interval
		d <- as.numeric(e) - as.numeric(s)

		# if min.coverage is not NA, each var of 'x' is tested to 
		# see if it respect the min.coverage defined.
		#
		# The test has to be done on each var individually because
		# NA values must be removed to calculate the effective
		# coverage.
		#
		# the result is kept in a vector
		if( !is.na(mc) ) {
			valid <- sapply(x@data, w=w, m=mc*d,
				function(var, w, m)
				{
					return( sum(w[!is.na(var)]) >= m )
				} )
		} else {
			valid <- rep(TRUE, length(x))
		}

		# init of the structure which will contain the results
		res <- list()

		# FUN is apply with defined arguments over each valid var
		# (against the min.coverage constraint)
		res[which(valid)] <- sapply(x@data[valid],
			      function(x, FUN, f.as)
			      {
				      f.as <- c(list(x), f.as)
				      do.call(FUN, f.as)
			      },
			      FUN, f.as)

		# var that are not valid (against the min.coverage constraint)
		# are defined and set to NA
		res[which(!valid)] <- NA

		res <- as.list(unlist( res ))

		# formatting the result to a TimeIntervalDataFrame
		names(res) <- names(x)
		res <- as.data.frame(res)
		res <- TimeIntervalDataFrame(s, e, data=res)
		return (res)
	}

	# apply of FUN and arguments over each splitted datas using 
	# the calc function defined above. 
	result <- mapply(calc, splitted, s.INDEX, e.INDEX,
			 MoreArgs=list(mc=min.coverage, ws=weights.arg,
				       FUN=FUN, f.as=fun.args),
			 SIMPLIFY=FALSE)

	# keep.INDEX indicates if initial var defined in INDEX must be 
	# kept in the final resulting list

	if( keep.INDEX )
		for( n in names(INDEX) )
		{
			# if var in INDEX and X have the same name,
			# those of INDEX are renamed
			new.n <- n
			i <- 1
			while( new.n %in% names(X) ) {
				new.n <- paste(n, i, sep='.')
				i <- i+1
			}
			if( new.n != n)
				warning(sprintf("'%s' in INDEX renamed as '%s'",
						n, new.n))

			# adding INDEX var in the result (eventually renamed)
			result <- mapply(function(x, i, value)
			       {
				       x@data[[i]] <- rep(value, nrow(x))
				       return( x )
	       		       },
			       x=result, value=INDEX[[n]],
			       MoreArgs=list(i=new.n), SIMPLIFY=FALSE)
		}

	# should the result be a list of TimeIntervalDataFrame or 
	# a unique TimeIntervalDataFrame ?

	if( simplify )
	{
		data <- sapply(result, function(x) unlist(as.data.frame(x)))
		if(is.null(dim( data )))
			data <- matrix(data, nrow=1,
				       dimnames=list(names(result[[1]])))

		data <- as.data.frame(t(data))

		start <- sapply(result, start)
		start <- as.POSIXct(start, origin=origin)
		end <- sapply(result, end)
		end <- as.POSIXct(end, origin=origin)

		result <- TimeIntervalDataFrame(start, end, data=data,
						timezone=timezone(INDEX))
	}
	
	return( result )
	} )

# tapply method to apply FUN over TimeIntervalDataFrame splitted by 
# a time period (POSIXctp) object.
# this method transform INDEX to a TimeIntervalDataFrame and then call the 
# main tapply method (see above).

setMethod (
	'tapply',
	signature(X='TimeIntervalDataFrame', INDEX='POSIXctp'),
	definition=function (X, INDEX, FUN, ..., min.coverage=1,
			     weights.arg=NULL, merge.X=TRUE,
			     split.X=FALSE, simplify=TRUE)
{

	s <- min(start(X))
	e <- max(end(X))
	tzone <- timezone( X )

	INDEX <- TimeIntervalDataFrame(s, e, period=INDEX, timezone=tzone)

	tapply( X, INDEX, FUN, ..., min.coverage=min.coverage, 
	       weights.arg=weights.arg, merge.X=merge.X, split.X=split.X,
	       simplify=simplify)
} )

# tapply method to apply FUN over TimeIntervalDataFrame splitted by 
# time intervals (POSIXcti) object.
# this method transform INDEX to a TimeIntervalDataFrame and then call the 
# main tapply method (see above).

setMethod (
	'tapply',
	signature(X='TimeIntervalDataFrame', INDEX='POSIXcti'),
	definition=function (X, INDEX, FUN, ..., min.coverage=1,
			     weights.arg=NULL, merge.X=TRUE,
			     split.X=FALSE, simplify=TRUE)
{
	tzone <- attributes(start(INDEX))$tzone[1]
	INDEX <- TimeIntervalDataFrame(start(INDEX), end(INDEX), tzone)

	tapply(X, INDEX, FUN, ..., min.coverage=min.coverage,
	       weights.arg=weights.arg, merge.X=merge.X,
	       split.X=split.X, keep.INDEX=FALSE, simplify=simplify)
} )

