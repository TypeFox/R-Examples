Meth <-
function( data = NULL,
          meth = "meth", item="item", repl=NULL, y="y",
         print = !is.null(data),
     keep.vars  =!is.null(data) )
{
dfr <- deparse(substitute(data))
was.dfr <- is.data.frame( data )
if( was.dfr ) dfr.nam <- names( data )

# Are measurements supplied in the wide format?
methods.in.y <- ( if( is.data.frame(y) ) ncol(y)>1 else FALSE  ) |
                (      is.character(y) & length(y)>1           ) |
                (        is.numeric(y) & length(y)>1 & was.dfr )

# Select the correct columns from a supplied data frame
if( was.dfr )
  {
  # Vector to collect the rows used constructing the Mehh object
  taken <- numeric(0)
  # Row numbers to use
  Nr <- nrow(data)
  rows <- 1:Nr

  ##############################################################################
  # First the case where several columns of the data argument constitute the
  # different methods
if( methods.in.y )
  {
  # Which columns contain the different methods
  if( is.character(y) )
    y <- match( y, dfr.nam )
  if( is.numeric(y) )
    {
    # Method names and the factor
    Mn <- names(data)[y]
    Nm <- length( Mn )
    meth <- rep( Mn, each=Nr )
    # Where were the y's
    taken <- y.col <- y
    # Put them under each other by the column-major default for matrices
    y <- as.vector( as.matrix(data[,y]) )
    }
  else stop( "y must be character or numeric if the data argument is given.\n" )

  # Items
  if( !missing(item) & is.character(item) )
    item <- match( item, dfr.nam )
  if( is.numeric(item) & length(item)==1 )
    {
    taken <- c(taken,item.col<-item)
    item <- data[,item]
    }
  else item <- rows

  # Replicates
  if( !missing(repl) & is.character(repl) )
    repl <- match( repl, dfr.nam )
  if( is.numeric(repl) & length(repl)==1 )
    {
    repl <- data[,repl]
    taken <- c(taken,repl.col<-repl)
    }
  else repl <- make.repl( data.frame(meth=rep(1,Nr),item=item) )$repl

  # Expand items and replicates
  if( length(item)==Nr ) item <- rep( item, Nm )
  if( length(repl)==Nr ) repl <- rep( repl, Nm )

  # To select the correct rows of other variables just transported across
  rows <- rep( rows, Nm )
  # End of several y-columns
  }

else
  {
  ##############################################################################
  # Then the case where we just select columns from the data argument

  # Method
  if( is.character(meth) )
    meth <- match( meth, dfr.nam )
  if( is.numeric(meth) & length(meth)==1 )
    {
    taken <- meth.col <- meth
    meth <- data[,meth]
    }
  if( is.na(meth)[1] ) stop( "\nmeth not properly specified.")

  # Item
  if( is.character(item) )
    item <- match( item, dfr.nam )
  if( is.numeric(item) & length(item)==1 )
    {
    taken <- c(taken,item.col<-item)
    item <- data[,item]
    }
  else item <- rows
  if( is.na(item)[1] ) stop( "\nitem not properly specified.")

  # Replicate
  if( is.null(repl) & "repl" %in% dfr.nam )
    repl <- "repl"
  if( is.character(repl) )
    repl <- match( repl, dfr.nam )
  if( is.numeric(repl) & length(repl)==1 )
    {
    taken <- c(taken,repl.col<-repl)
    repl <- data[,repl]
    }
  else repl <- rep(1,nrow(data))

  # Measurements
  if( is.character(y) )
    y <- match( y, dfr.nam )
  if( is.numeric(y) & length(y)==1 )
    {
    taken <- c(taken,y.col<-y)
    y <- data[,y]
    }
  if( is.na(y)[1] ) stop( "\ny not properly specified.")

  # End of one y-column
  }

# End of the data.frame case
}

# The following also covers the situation where y is supplied as a dataframe
if( length(y) > 1 &
    is.list(y) )
  {
  Mn <- if( !is.null(names(y)) ) names(y)
        else paste( "Method", 1:length(y), sep="" )
  Nm <- length( Mn )
  Nr <- length( y[[1]] )
  meth <- rep( Mn, each=Nr )
  if( missing(item) ) item <- 1:Nr
  item <- rep( item, Nm )
  if( missing(repl) ) repl <- rep(1,Nr)
  repl <- rep( repl, Nm )
     y <- unlist(y)
  }

# The resulting dataframe
res <- data.frame( meth = factor(meth),
                   item = factor(item),
                   repl = factor(repl),
                      y = y )

# Transfer other columns of the dataframe:
# Keep all other variables
if( was.dfr & is.logical(keep.vars) )
  {
  if( keep.vars )
    {
    res <- data.frame( res, data[rows,-taken] )
    names(res)[-(1:4)] <- dfr.nam[-taken]
    }
  }
# Keep only the selected variables
if( was.dfr & is.character(keep.vars) )
  keep.vars <- match( keep.vars, names(data) )
if( was.dfr & is.numeric(keep.vars) )
  {
  res <- data.frame( res, data[rows,keep.vars] )
  names(res)[-(1:4)] <- dfr.nam[keep.vars]
  }

# Remove missing y-values (and possible factor levels only in these)
if( any(is.na(res$y)) ) res <- res[!is.na(res$y),]

# Check if there actually are replicates even if not indicated
made.repl <- FALSE
if( any( table(res$meth,res$item)> 1 ) )
  if( length(table(res$repl))==1 )
    if( names(table(res$repl))=="1" )
      {
      res <- make.repl( res )
      made.repl <- TRUE
      }
class( res ) <- c("Meth","data.frame")

if( print & was.dfr )
cat( "The following variables from the dataframe\n\"",
     dfr, "\" are used as the Meth variables:",
     paste( if( exists("meth.col") ) paste("\nmeth:",dfr.nam[meth.col]),
            if( exists("item.col") ) paste("\nitem:",dfr.nam[item.col]),
            if( exists("repl.col") ) paste("\nrepl:",dfr.nam[repl.col]),
                                           "\n   y:",
                                               paste(dfr.nam[   y.col],collapse=" "),
            "\n" ),
     sep="" )

if( print )  print( summary.Meth( res ) )

if( made.repl )
  cat( "\nNOTE: Replication numbers generated in the order of the data\n" )
if( max( with( res, table( meth, item, repl ) ) ) > 1 )
  cat( "\nWARNING: You have chosen a replicate variable which is not unique\n",
       "        within each (meth,item).\n" )

attr(res,"row.names") <- 1:nrow(res)
invisible( res )
}

# Calculation of mean/median/min/max over replicates
mean.Meth <-
function( x, na.rm=TRUE, simplify=TRUE, ... )
{
tmp <- aggregate( x$y,
                  list( meth=x$meth, item=x$item ),
                  FUN=mean,
                  na.rm=na.rm, ... )
if( simplify )
  {
  names( tmp )[3] <- "y"
  tmp
  }
else
  {
  FUN.name <- deparse( substitute(FUN) )
  names( tmp )[3] <- paste( FUN.name, "y", sep="." )
  tmp <- merge(x,tmp)
  }
return( invisible( Meth( tmp, print=FALSE ) ) )
}

sort.Meth <-
function( x, ... )
{
tmp <- x[order(x$meth,x$item),]
tmp$sort.y <- ave( tmp$y, list(tmp$meth,tmp$item), FUN=sort )
tmp <- merge(x,tmp)
return( invisible( Meth( tmp, print=FALSE ) ) )
}

# Utilities needed to preserve the Meth attribute
subset.Meth <-
function( x, ... )
{
y <- base::subset.data.frame(x, ...)
# In order to reduce the factor levels to the ones actually present, use Meth()
return( Meth(y,keep.vars=TRUE,print=FALSE) )
}

transform.Meth <-
function( `_data`, ...)
{
    save.at <- attributes(`_data`)
    y <- base::transform.data.frame(`_data`, ...)
    save.at[["names"]] <- attr(y, "names")
    attributes(y) <- save.at
    y
}
