
setMethod( "[", signature( x = "groupVec", i = "ANY" ),
  function(x, i, ..., drop = TRUE )
  {
    # subscripting for groupVec
    # drop argument is ignored
    x@columns <- lapply( x@columns, "[", i, drop = FALSE )
    x
  })

setReplaceMethod( "[",
  signature( x = "groupVec", i="ANY", j="ANY", value = "groupVec" ),
  function( x, i, j, ..., value )
  {
    # subscript replace for groupVec when RHS is also a groupVec
    # calls the subscript replace with RHS a list
    value <- as( value, class(x))
    if( !identical( value@names, x@names ))
      stop( "Cannot replace -- different structures" )
    len = length(value@names)
    for(k in 1:len){
      x@columns[[k]][i] <- value@columns[[k]][...]
    }
    x
  })


setReplaceMethod( "[", signature( x = "groupVec",  i="ANY", j="ANY", value = "list" ),
   function( x, i, j, ..., value )
   {
  # subscript replace for groupVec when RHS is a list
     CurrData <- x@columns
     DataWidth <- length( CurrData )
  # verify that the value list is same length as data list
     if( DataWidth != length( value ))
       stop( "Replacement value has wrong structure -- cannot assign" )
     if( !DataWidth )
       return( x )
     ls <- sapply( value, "length" )
     l <- ls[1]
     if( !all( ls == l ))
      stop( "Replacement value columns have different lengths" )

     CurrData <-
       lapply( 1:DataWidth,
              function( n, x, value, cls, i )
              {
                ## apparently, x[[n]][...] is not allowed on the left side
                ## of an assignment -- all versions of S complain
                tmp <- x[[n]]
                tmp[i] <- value[[n]]
                if( !is( tmp, cls[n] ))
                  stop( "Replacement value has wrong class")
                tmp
              },
              CurrData, value, x@classes, i )
     x@columns <- CurrData
     x
   })


setReplaceMethod( "[", signature( x = "groupVec",  i="ANY", j="ANY", value = "ANY" ),
  function( x, i, j, ..., value )
  {
    # subscript replace for groupVec when RHS is anything but
    # groupVec or list
    # coerce the RHS to the class of the LHS, and subscript
    value <- as( value, class( x ))
    x[i] <- value
    x
  })

setMethod( "[[", signature( x = "groupVec" , i = "ANY" ),
  function(x, i, ... )
  {
    x[i]
  })

setReplaceMethod( "[[", signature( x = "groupVec" , i = "ANY"),
  function( x, i, ..., value )
  {
    x[i] <- value
    x
  })

setMethod( "length<-", signature( x = "groupVec" ),
         function( x, value )
       {
         # length replacement for groupVec
         x@columns <- lapply( x@columns, "length<-", value )
         x
       })


setMethod( "length", signature( x = "groupVec" ),
       function( x )
       {
       # length for groupVec
	 mydata <- x@columns
	 if( length( mydata ) < 1 )
	   return( 0 )
	 length( mydata[[1]] )
       })

setMethod( "c", signature(x="groupVec"),
  function(x, ...){
    ## Concatenates the columns slot of the two groupVec objects, if
    ## they are of the same class or can be coerced to be the same.
    ## can be called with one or the other of x, y a groupVec class --
    ## otherwise it is likely to fail.
    arglist <- list(...)
    if(length(arglist)==0) return(x)
    lens <- sapply(arglist, length)
    if(!any(lens > 0)) {
      return(x)
    }
    arglist <- arglist[lens > 0]
    if(length(arglist) > 1)
      c(c(x, arglist[[1]]), do.call("c", arglist[-1]))
    else {
      if(!identical(names(arglist[[1]]), names(x)))
        stop("cannot concatenate -- all slots must have the same names")
      if(!identical(arglist[[1]]@classes, x@classes))
        stop("cannot concatenate -- all slots must have the same classes")
      y <- as( arglist[[1]], class( x ))
      l = length(x@columns)
      for(n in 1:l){
        yn = as(y@columns[[n]], class(x@columns[[n]]))
        x@columns[[n]] <- c(x@columns[[n]], yn)
      }
      x
    } 
  })

setMethod( "is.na", "groupVec",
     function( x )
     {
      # is.na for groupVec
      # returns a logical vector whose elements correspond to the rows of x,
      # with each value having the value any( is.na( columns of x ))
       dat <- x@columns
       len <- length( dat )
       if(( len < 1 ) || ( length( x ) < 1 ))
	 return( logical( 0 ))
       res <- FALSE
       for( i in 1:len )
	 res <- res | is.na(dat[[i]])
       res
     })

setMethod( "is.infinite", "groupVec",
     function( x )
     {
      # is.infinite for groupVec
      # returns a logical vector whose elements correspond to the rows of x,
      # with each value having the value any( is.infinite( columns of x ))
       dat <- x@columns
       len <- length( dat )
       if(( len < 1 ) || ( length( x ) < 1 ))
	 return( logical( 0 ))
       res <- FALSE
       for( i in 1:len )
	 res <- res | is.infinite(dat[[i]])
       res
     })

setMethod( "is.finite", "groupVec",
     function( x )
     {
      # is.finite for groupVec
      # returns a logical vector whose elements correspond to the rows of x,
      # with each value having the value any( is.finite( columns of x ))
       dat <- x@columns
       len <- length( dat )
       if(( len < 1 ) || ( length( x ) < 1 ))
	 return( logical( 0 ))
       res <- FALSE
       for( i in 1:len )
	 res <- res | is.finite(dat[[i]])
       res
     })

setMethod( "is.nan", "groupVec",
     function( x )
     {
      # is.nan for groupVec
      # returns a logical vector whose elements correspond to the rows of x,
      # with each value having the value any( is.nan( columns of x ))
       dat <- x@columns
       len <- length( dat )
       if(( len < 1 ) || ( length( x ) < 1 ))
	 return( logical( 0 ))
       res <- FALSE
       for( i in 1:len )
	 res <- res | is.nan(dat[[i]])
       res
     })

setMethod( "unique", signature( x = "groupVec" ),
          function( x, ... )
          {
            dups = duplicated(x)
            if(any(dups)){
              x@columns = lapply(x@columns,
                function(col, drop) col[!drop], dups)
            }
            return(x)
          })

setMethod( "duplicated", "groupVec",
     function( x, incomparables=FALSE )
     {
      ## duplicated for groupVec
      ## returns a logical vector whose elements correspond to the rows of x,
      ## with each value having the value all( duplicated( columns of x ))
       dat <- x@columns
       len <- length( dat )

       if(( len < 1 ) || ( length( x ) < 1 ))
	 return( logical( 0 ))

       # now find the unique items in each column
       uniques <- lapply( dat, "unique" )
       # and match them in each column to make categories
       cats <- lapply( 1:len, function( col, data, uniq )
		      match( data[[col]], uniq[[col]] ), dat, uniques )
       # now we have each column reduced to integer categories
       # paste each row together and use duplicated on the
       # result
       ret <- duplicated( do.call( "paste", cats ))

       # and take care of the incomparables argument
       if( !identical( incomparables, FALSE ))
       {
	 whereinc <- ( match( x, incomparables, nomatch=-1) > 0 )
	 ret[whereinc] <- FALSE
       }
       ret
     })

setMethod( "show", "groupVec", function( object )
{
  # show function for groupVec
  mynames <- names( object )
  myclasses <- object@classes
  objSlots <- names( getSlots( class(object) ))
  GVslots <- names( getSlots( "groupVec" ))
  objSlots <- objSlots[ is.na( match( objSlots, GVslots )) ]

  cat( "An object of class \"", class( object ), "\"\n", sep = "" )
  if( length( mynames ))
    for( i in 1:length( mynames ))
    {
      cat( "column \"", mynames[i], "\" (", myclasses[i], "):\n", sep = "" )
      show( groupVecColumn( object, mynames[i] ))
    }
  else
    cat( "no columns, no data\n" )

  for( slname in objSlots )
  {
    cat( "slot \"", slname, "\":\n", sep = "" )
    show( slot( object, slname ))
  }
})

setMethod( "summary", "groupVec", function( object, ... )
{
  # summary function for groupVec
  mynames <- names( object )
  myclasses <- object@classes
  allslots <- getSlots( class(object) )

  GVslots <- names( getSlots( "groupVec" ))
  whereslots <- is.na( match( names( allslots ), GVslots ))
  slotnames <- names(allslots)[whereslots]
  slotclasses <- allslots[whereslots]

  sumry <- array( "", c( length( slotnames ) + length( mynames ), 3 ),
		  list( c( mynames, slotnames ),
		        c( "Slot/Column", "Length", "Class" )))

  if( length( mynames ))
  {
    sumry[mynames,1] <- "Column"
    sumry[mynames,2] <- format( length( object ))
    sumry[mynames,3] <- myclasses
  }

  if( length( slotnames ))
  {
    sumry[slotnames,1] <- "Slot"
    sumry[slotnames,3] <- slotclasses

    for( slname in slotnames )
      sumry[slname,2] <- format( length(slot(object,slname)))
  }

  oldClass(sumry) <- "table"
  sumry
})

setMethod( "all.equal.list", "groupVec",
function (target, current, ...) 
{
    check.attributes <- list(...)[["check.attributes"]]
    if (is.null(check.attributes)) check.attributes <- TRUE
    msg <- if (check.attributes) 
        attr.all.equal(target, current, ...)
    iseq <- if ((length(target) == length(current)) &&
        all(match(names(target), names(current), nomatch=0) > 0)) {
          names(target)
        } else {
          if (!is.null(msg)) {
            msg <- msg[-grep("\\bNames\\b", msg)]
            msg <- c(msg, paste("Name mismatch: comparison on intersection of names"))
          } 
          intersect(names(target), names(current))
        }
    for (i in iseq) {
        mi <- all.equal(groupVecColumn(target,i), groupVecColumn(current,i), 
            ...)
        if (is.character(mi)) 
            msg <- c(msg, paste("Component ", i, ": ", mi, sep = ""))
    }
    if (is.null(msg)) 
        TRUE
    else msg
})

setMethod("names", "groupVec", function(x){
  x@names
})

setMethod("rep", signature(x = "groupVec"),
function(x, ...)
{
  newCols = lapply(x@columns, function(y) rep(y, ...))
  x@columns = newCols
  x
})
