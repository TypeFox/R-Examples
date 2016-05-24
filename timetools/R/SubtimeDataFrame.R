# definition de la classe
#------------------------
setClass (Class = 'SubtimeDataFrame', 
	  representation = representation (when='POSIXst',
					   data='data.frame'),
	  prototype = prototype (when=new('POSIXst'),
				 data=data.frame()),
	  validity=function(object) {
		  if (length (when (object)) != nrow (object))
			  stop ("In a 'SubtimeDataFrame, 'data' must have a number of rows as long as 'when'.")
		  return (TRUE)
	  } )

# constructeurs
#--------------
SubtimeDataFrame <- function (when, data=NULL, ...) {
	if (is.null (data))
		data <- data.frame (matrix (NA, ncol=0, nrow=length(when) ) )
	new ('SubtimeDataFrame', when=when, data=data)
}

# definition des accesseurs de l'objet
#-------------------------------------
setMethod (f='when', signature='SubtimeDataFrame',
	   definition=function(x, ...) return(x@when) )

setMethod (f='timezone', signature='SubtimeDataFrame',
	   definition=function(object) return(timezone(when(object))) )

unit.SubtimeDataFrame <- function(x, ...) unit(when(x))

of.SubtimeDataFrame <- function(x, ...) of(when(x))

# mise en forme pour / et affichage
#----------------------------------
print.SubtimeDataFrame <- function (x, ...)
	print(data.frame (when=format(when(x), ...), x@data) )

setMethod ('show', 'SubtimeDataFrame',
	   function (object) {
	   print(data.frame (when=format(when(object)), object@data) )
	   } )
tail.SubtimeDataFrame <- function (x, ...)
		print (tail (data.frame (when=format(when(x), ...), x@data) ) )

head.SubtimeDataFrame <- function (x, ...)
		print(head (data.frame (when=format(when(x), ...), x@data) ) )

summary.SubtimeDataFrame <- function (object, ...)
		print (summary (data.frame (when=format(when(object), ...),
					    object@data) ) )
# format

# defintion des accesseurs aux donnees
#-------------------------------------
'[.SubtimeDataFrame' <- function(x, i, j, drop=FALSE) {
	n.args <- nargs() - hasArg(drop)
	if (missing (j) & n.args==2) {
		j <- i
		i <- seq_len(nrow(x))
	}
	if(missing(i)) i <- seq_len(nrow(x))
	y <- new ('SubtimeDataFrame', 
	     when =when (x)[i],
	     data = x@data[i, j, drop=drop])
	validObject(y)
	return(y)
}
setMethod (f='[[', signature='SubtimeDataFrame',
	   definition=function(x, i, ...) {
		   '[[.data.frame'(x@data, i, ...)
	   })
setMethod (f='$', signature='SubtimeDataFrame',
	   definition=function(x, name) {
		   do.call ('$', list(x=x@data, name=name))
	   })

'[<-.SubtimeDataFrame' <- function(x, i, j, value) {
	n.args <- nargs()
	if (missing (j) & n.args==3) {
		j <- i
		i <- seq_len(nrow(x))
	}
	if(missing(i)) i <- seq_len(nrow(x))
	x@data[i,j] <- value
	validObject(x)
	return(x)
}
'[[<-.SubtimeDataFrame' <- function(x, i, j, value) {
   if (missing (j) )
	   x@data[[i]] <- value else
	   x@data[[i,j]] <- value
   validObject(x)
   return(x)
}
setMethod (f='$<-', signature='SubtimeDataFrame',
	   definition=function(x, name, value) {
		   x@data <- "$<-.data.frame"(x@data, name, value)
		   validObject(x)
		   return(x)
	   })

setMethod (f='dim', signature='SubtimeDataFrame',
	   definition=function(x) dim (x@data))
setMethod (f='nrow', signature='SubtimeDataFrame',
	   definition=function(x) nrow (x@data))
setMethod (f='ncol', signature='SubtimeDataFrame',
	   definition=function(x) ncol (x@data))
row.names.SubtimeDataFrame <- function(x) row.names (x@data)
'row.names<-.SubtimeDataFrame' <- function(x, value)
{
	row.names (x@data) <- value
	x
}
setMethod (f='names', signature='SubtimeDataFrame',
	   definition=function(x) names (x@data))
setMethod (f='names<-', signature='SubtimeDataFrame',
	   definition=function(x, value) {
		   names (x@data) <- value
		   x
	   } )

# Math

# manipulation
#-------------
# # fonction réalisée en S3 pour ne pas imposer de 'signature'
# rbind.SubtimeDataFrame <- function (...) {
# 	dots <- list (...)
# 	names(dots) <- NULL
# 	dots <- lapply (dots, function(x) as(x, 'SubtimeDataFrame') )
# 	when <- unlist (lapply (dots, when) )
# 	df <- do.call("rbind", lapply(dots, function(x) x@data) )
# 	new('SubtimeDataFrame', when=when, data=df)
# }
# cbind # a faire eventuellement entre un Time*DataFrame et une data.frame
merge.SubtimeDataFrame <- function(x, y, by, all=TRUE, ...) {
		if (!inherits(y, 'SubtimeDataFrame'))
			stop ("'y' must be a 'SubtimeDataFrame'.")
		if( any(duplicated(when(x))) | any(duplicated(when(y))) )
			stop("'when' slots must be unique in each SubtimeDataFrame")
		if (missing (by) ) by <- NULL

		when.vec <- list (when(x), when(y))

		if( unit(x) != unit(y) )
			stop("x and y must have same unit") else
			u <- unit(x)
		if( of(x) != of(y) )
			stop("x and y must have same of") else
			o <- of(x)
		if( timezone(x) != timezone(y) )
			stop("x and y must have same timezone") else
			tz <- timezone(x)

		x.data <- data.frame(when=as.numeric(format(when(x), '%v')),
				     x@data)
		y.data <- data.frame(when=as.numeric(format(when(y), '%v')),
				     y@data)
		z <- merge (x.data, y.data,
			    by=unique (c('when', by) ), all=all, ...)
		
		when <- POSIXst( z$when, u, o, tz )

		z <- new ('SubtimeDataFrame',
			  when=when,
			  data=z[setdiff(names(z), c('when'))])
		return (z)
	   }

setMethod ('lapply', signature('SubtimeDataFrame', 'ANY'),
	   function (X, FUN, ...)
	   {
		   res <- lapply (data.frame(X), FUN, ...)
		   if (all (sapply (res, length) == nrow(X))) {
			   X@data <- data.frame (res[names(X)])
		   } else {
			   warning ("Result has a number of rows differents from object. A data.frame is returned.")
			   X <- data.frame (res)
		   }
		   return (X)
	   } )


# acces/modification de certaines propriétés
#-------------------------------------------


# transformateur de classe
#-------------------------
# as.TimeInstantDataFrame.TimeIntervalDataFrame <- function(from, how=c('mid', 'start', 'end'), ...) {
#         how <- match.arg (how)
#         instant <- switch (how,
#                            mid=.POSIXct(rowMeans(data.frame(unclass(start(from)), unclass(end(from)))), tz=from@timezone),
#                            start=start(from), end=end(from))
#         to <- new ('TimeInstantDataFrame', instant=instant, data=from@data)
#         validObject(to)
#         return (to)
# }

as.data.frame.SubtimeDataFrame <- function (x, row.names=NULL, optional=FALSE,
			include.dates=FALSE, ...) {
	if (include.dates)
		return (data.frame (when=when(x), x@data) ) else
		return (x@data)
}
