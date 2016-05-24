# definition of the generic (S4) function
setGeneric (name='split')

# the following (S3) methods are wrappers to the split.default method
# defined in the package 'base'.
# Those methods allow to split POSIXctp, POSIXst  POSXcti objects
# which are similar to vectors.

split.POSIXctp <- function(x, f, drop=FALSE, ...)
{
	i <- seq_len(length(x))
	i <- split(i, f, drop)
	lapply(i, function(i, x) x[i], x)
}

split.POSIXst <- function(x, f, drop=FALSE, ...)
{
	i <- seq_len(length(x))
	i <- split(i, f, drop)
	lapply(i, function(i, x) x[i], x)
}

split.POSIXcti <- function(x, f, drop=FALSE, ...)
{
	i <- seq_len(length(x))
	i <- split(i, f, drop=drop)
	lapply(i, function(i, x) x[i], x)
}

# the following (S3) methods are wrappers to the split.data.frame
# method defined in the package 'base'.
# Those methods allow to split SubtimeDataFrame, TimeInstantDataFrame
# and TimeIntervalDataFrame in the same way that data.frame are
# used to be splitted.

split.SubtimeDataFrame <- function(x, f, drop=FALSE, ...)
{
	vect <- seq_len(nrow(x))
	w <- split (when(x), f, drop)
	data <- split (x@data, f, drop)
	x <- mapply (SIMPLIFY=FALSE, new, 'SubtimeDataFrame',
		     when=w, data=data, USE.NAMES=FALSE)
	names( x ) <- names( data )
	x
}

split.TimeInstantDataFrame <- function(x, f, drop=FALSE, ...)
{
	vect <- seq_len(nrow(x))
	i <- split (when(x), f, drop)
	data <- split (x@data, f, drop)
	x <- mapply (SIMPLIFY=FALSE, new, 'TimeInstantDataFrame',
		     instant=i, data=data, timezone=timezone(x),
		     USE.NAMES=FALSE)
	names( x ) <- names( data )
	x
}

split.TimeIntervalDataFrame <- function(x, f, drop=FALSE, ...)
{
	vect <- seq_len(nrow(x))
	s <- split (start(x), f, drop)
	e <- split (end(x), f, drop)
	data <- split (x@data, f, drop)
	x <- mapply( SIMPLIFY=FALSE, new, 'TimeIntervalDataFrame',
		     start=s, end=e, data=data, timezone=x@timezone,
		     USE.NAMES=FALSE)
	names( x ) <- names( data )
	x
}

# Since POSIXst, POSIXctp and POSIXcti objects are simila to vector
# it must be possible to split other type of objects against those ones.
# Here are the definition of such methods

setMethod('split', signature('ANY', 'POSIXst'),
	  function(x, f, drop=FALSE, ...)
	  {
		  f <- as.numeric(f)
		  split(x, f, drop=FALSE, ...)
	  } )

setMethod('split', signature('ANY', 'POSIXctp'),
	  function(x, f, drop=FALSE, ...)
	  {
		  f <- format(f)
		  split(x, f, drop=FALSE, ...)
	  } )

setMethod('split', signature('ANY', 'POSIXcti'),
	  function(x, f, drop=FALSE, ...)
	  {
		  f <- format(f, ...)
		  split(x, f, drop=FALSE)
	  } )


# the next methods are more specialised methods and each one will
# be described independently.

# split a TimeIntervalDataFrame into another TimeIntervalDataFrame.
# The method take each time interval of the first TimeIntervalDataFrame
# (TitDF) and search in which time intervals of the second it
# intersects. Each time interval of the first TItDF can intersect with
# none, one or several time intervals of the second TItDF. The arguments
# 'split.x' is defined to tell the method what to do : 
# - if the time interval in the first TItDF (ti1) doesn't match any in the 
# second TItDF, nothing to do ;
# - if it (ti1) matches one in the second TItDF (ti2) and is
# included inside it, it (ti1) is entirely taken in the final result ;
# - if it (ti1) intersects one and only one (ti2) inside the second
# TItDF, (ti1) is truncated to be included inside (ti2) if 'split.x' is
# TRUE and (ti1) is removed if 'split.x' is FALSE ;
# - if it (ti1) is over several time intervals of the second TItDF
# (ti2.a, ti2.b, etc.) :
# -- if 'split.x' is TRUE, (ti1) is truncated into each ti2.x to be
# included inside each one ;
# -- if 'split.x' is FALSE, (ti1) is removed.

setMethod('split',
	  signature=signature(x='TimeIntervalDataFrame',
			      f='TimeIntervalDataFrame'),
	  definition=function(x, f, ...,
			      split.x=FALSE, keep.f=TRUE)
{
	if (timezone(x) != timezone (f))
	{
		warning("'x' and 'f' have a different timezone. The timezone of 'f' is taken for the result")
		timezone(x) <- timezone(f)
	}

	#=========================================================================
	if (homogeneous (x) && continuous(x) &&
	    homogeneous (f) && continuous(f) &&
	    (as.numeric(period(f)) / as.numeric(period(x)))%%1 == 0 &&
	    start(x)[1] == start(f)[1] ) {
	# trivial case where time support of 'x' is a subdivision of the time
	# support of 'f', a specific algorithm is used to improve speed
	# calculation
	message ("Simplified splitting algorithm")

	# x and f are ordered and the number of rows of x per rows of f
	# is calculated
	x <- x[order(start(x)), ]
	f <- f[order(start(f)), ]

	nb.x <- (as.numeric(period(f)) / as.numeric(period(x)))

	# definition of a vector containing the grouping value
	split.vect <- rep(1:nrow(f), each=nb.x)

	# every needed data to build the final result are splitted over
	# the vector defined above
	s <- split(start(x), split.vect)
	e <- split(end(x), split.vect)
	tz <- timezone(x)
	data <- split(x@data, split.vect)

	# building of the result
	result <- mapply(
		function(s, e, d, tz)
			TimeIntervalDataFrame(start=s, end=e, timezone=tz, data=d),
    		s, e, data, MoreArgs=list(tz=tz))

	#=========================================================================
	} else {
	# universal algorithm which works for every case (for trivial case 
	# above to, but less efficiently)

	# retrieving informations from input data to define the arguments
	# to pass to the C functions (see below)

	int.x <- when (x)
	int.f <- when (f)
	   
	s.x <- start (int.x)
	e.x <- end (int.x)
	s.f <- start (int.f)
	e.f <- end (int.f)

	# the underlying C function determine the number of time
	# intersections between time intervals of 'x' and time intervals
	# of 'f'. This number is necessary for the enxt step.

	nb <- .C ('project_nb_intersections',
			as.integer(s.x), as.integer(e.x),
			as.integer(length(s.x)),
			as.integer(s.f), as.integer(e.f),
			as.integer(length(s.f)),
			nb=integer(1),
			NAOK=FALSE, PACKAGE='timetools')$nb

	# The 'whiches'  structure 
	# is a data.frame containing 3 variables. Each row of the
	# data.frame represents a time intersection between 'x'
	# and 'f'.
	# pos.x: the row index in 'x' corresponding to the time
	#	interserction ;
	# pos.f: the row index in 'f' corresponding to the time
	#	interserction ;
	# weight: this one is not used is this function.

	if (nb > 0) {

		whiches <- .C ('project_pos_weight',
			as.integer(s.x), as.integer(e.x),
			as.integer(length(s.x)),
			as.integer(s.f), as.integer(e.f),
			as.integer(length(s.f)),
			pos.x=integer(nb), pos.f=integer(nb),
			weight=integer(nb),
			NAOK=FALSE, PACKAGE='timetools')[c('pos.x', 'pos.f')]
		whiches <- as.data.frame (whiches)

	} else {

		whiches <- data.frame (pos.x=numeric(), pos.f=numeric())

	}

	# if 'split.x' is FALSE, time intervals of 'x' are not allowed to
	# be over several time intervals of 'f' neither to intersect one
	# time interval of 'f' being not completely included in it.
	# The first condition implies that in the 'pos.x' variable of 
	# 'whiches', every value must appear only once. Rows in 'whiches'
	# where the 'pos.x' value is not unique are removed.
	# WARNING : the second condition is not yet tested so we have
	# to keep it in mind to treat it later.
	#
	# 'splitted', 's' and 'e' (start and end of the time interval
	# of 'x') ans 'si' and 'ei' (start and end of the time interval
	# of 'f') variables are added to 'whiches'.

	if (!split.x) {
		whiches <- whiches[!whiches$pos.x %in%
			unique( whiches$pos.x[duplicated(whiches$pos.x)] ),]
		whiches$splitted <- rep(FALSE, nrow(whiches))
	} else {
		whiches$splitted <- whiches$pos.x %in%
			unique( whiches$pos.x[duplicated(whiches$pos.x)] )
	}

	whiches$s <- as.numeric(start(x))[whiches$pos.x]
	whiches$e <- as.numeric( end (x))[whiches$pos.x]

	whiches$si <- as.numeric(start(f))[whiches$pos.f]
	whiches$ei <- as.numeric( end (f))[whiches$pos.f]

	if (!split.x) {
		# with the new variables, we can now easily test the
		# second condition about time interval of 'x' that are
		# not commpletely included in a time interval of 'f'.
		# The condition is not met if the "x"'s interval begin before
		# or end after the "f"'s one.

		whiches <- whiches[whiches$s >= whiches$si &
				   whiches$e <= whiches$ei,]

	} else {
		# as explained in the description of the function,
		# if splitting of "x"'s intervals is allowed, those
		# that are not completely included in "f"'s intervals
		# are truncated : there start ('s') and end ('e')
		# are modified in accordance with this rule.

		early  <- whiches$s < whiches$si
		lately <- whiches$e > whiches$ei
		whiches$s[early ] <- whiches$si[ early]
		whiches$e[lately] <- whiches$ei[lately]
	}

	# some of the time intervals of 'f' may have intersection with
	# none of the time intervals of 'x'.
	# First, we determine the TimeIntervalDataFrame for those
	# that are not empty.

	# get 'pos.f' values for not empty "f"'s intervals and
	# corresponding start and end.
	f.notempty <- sort(unique( whiches$pos.f ))
	start <- as.numeric(s.f[f.notempty])
	end <- as.numeric(e.f[f.notempty])

	# get whiches splitted into "f"'s intervals
	whiches <- split (whiches, whiches$pos.f)
 
	# split effectif des donnees
	# pour chaque element du whiches (et donc de start et end)
	# on créé un TimeIntervalDataFrame

	tz <- timezone(x)

	result <- list()
	result[f.notempty] <- mapply(
		function(w, s, e, x, tz) {

			sx <- w$s
			ex <- w$e
			if( any(w$splitted) ) {
				sx[w$splitted] <- sapply(sx[w$splitted], max, s)
				ex[w$splitted] <- sapply(ex[w$splitted], min, e)
			}

			TimeIntervalDataFrame(start=as.POSIXct(sx, origin=origin),
					      end=as.POSIXct(ex, origin=origin),
					      timezone=tz,
			      		      data=x@data[w$pos.x,,drop=FALSE])
		},
		whiches, start, end,
		MoreArgs=list(x, tz=tz), SIMPLIFY=FALSE)

	# for empty "f"'s intervals, empty TimeIntervalDataFrame are
	# set in the resultint structure.

	result[setdiff(1:nrow(f), f.notempty)] <- TimeIntervalDataFrame(
		as.POSIXct(character()),
		as.POSIXct(character()),
		data = x@data[0,,drop=FALSE])

	} # end of switching over trivial case or not
	#=========================================================================

	# adding initial values of 'f' in the result if asked.

	if( keep.f )
		for( n in names(f) )
		{
			# if 'x' and 'f' have variables with identical names
			# those of 'f' are modified
			new.n <- n
			i <- 1
			while( new.n %in% names(x) ) {
				new.n <- paste(n, i, sep='.')
				i <- i+1
			}
			if( new.n != n)
				warning(sprintf("'%s' in f renamed as '%s'",
						n, new.n))

			# add "f"'s vars
			result <- mapply(function(x, i, value)
			       {
				       x@data[[i]] <- rep(value, nrow(x))
				       return( x )
	       		       },
			       x=result, value=f[[n]],
			       MoreArgs=list(i=new.n), SIMPLIFY=FALSE)
		}

	return( result )
} )

# split a TimeIntervalDataFrame into a time period (of length 1).
# A TimeIntervalDataFrame is created (cf TimeIntervalDataFrame
# constructor) and the the above method is called.

setMethod (
	'split',
	signature(x='TimeIntervalDataFrame', f='POSIXctp'),
	definition=function (x, f, ...,  split.x=FALSE)
{
	s <- min(start(x))
	e <- max(end(x))
	tzone <- timezone( x )
	f <- TimeIntervalDataFrame(s, e, period=f, timezone=tzone)

	split(x, f, ..., split.x=split.x)
} )

# split a TimeIntervalDataFrame into time intervals (POSIXcti).
# It is exactly the same as splitting a TimeIntervalDataFrame into
# another except that 'f' has not data.
# So a TimeIntervalDataFrame is created according to 'f' and the 
# the above method is called.

setMethod('split',
	  signature=signature(x='TimeIntervalDataFrame', f='POSIXcti'),
  	  definition=function(x, f, ..., split.x=FALSE)
{
	tzone <- attributes(start(f))$tzone[1]
	f <- TimeIntervalDataFrame(start(f), end(f), tzone)

	split(x, f, ..., split.x=split.x)
} )

# TODO : split,TimeInstantDataFrame,TimeIntervalDataFrame-method

