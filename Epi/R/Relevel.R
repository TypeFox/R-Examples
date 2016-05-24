
# The levels method is already defined (in the utils package)
# and hence imported in the NAMESPACE file
levels.Lexis <-
function( x )
{
union( levels(x$lex.Cst), levels(x$lex.Xst) )
}

# The Relevel method
Relevel <- function (x, ...) UseMethod("Relevel")

# The factor method is the default method
Relevel.default <-
Relevel.factor <-
  function( x, ref, first=TRUE, collapse="+", ... )
  {
  # Function that collapses multiple sets of levels of a factor
  #
  if( !is.factor(x) )
    {
    argnam <- deparse( substitute(x) )
    f <- factor( x )
    cat( "WARNING: ", argnam,
         "has been converted to a factor with levels:\n",
         levels( f ) )
    }
  else
    f <- x

  # This is a copy of the relevel function from the base package:
  #
  relev <- function (x, ref, ...)
  {
    lev <- levels(x)
    if ( is.character( ref ) )
         ref <- match(ref, lev)
    if ( any( is.na( ref ) ) )
         stop( "any values in ref must be an existing level" )
    nlev <- length( lev )
    if ( any( ref < 1 ) || any( ref > nlev ) )
         stop( paste( "ref=", paste( ref, collapse="," ),
                      ": All elements must be in 1:", nlev, sep="" ) )
    factor(x, levels = lev[c(ref, seq(along = lev)[-ref])])
  }

  # If called with a non-list argument assume reshuffling of levels
  #
  if( !is.list( ref ) )
    fnew <- relev( f, ref )

  # If called with a list collapse levels in each list element.
  #
  if( is.list( ref ) )
    {
      fnew <- f
      newnames <- levels( f )
      uninames <- character( length( ref ) )
      for( s in 1:length(ref) )
         if ( is.character( ref[[s]] ) ) ref[[s]] <- match( ref[[s]], levels(f) )
      # Check for replicates in levels to be grouped
      if( any( (tt<-table(unlist(ref))) > 1 ) )
        stop("Factor level(s) ", levels( f )[as.numeric(names(tt)[tt>1])],
             " has been allocated to more than one new level.")
      for( s in 1:length(ref) )
         {
         uninames[s] <- if( is.null( names( ref ) ) )
                          {
                          paste( levels( f )[ref[[s]]], collapse=collapse )
                          }
                        else if( names( ref )[s]=="" )
                          {
                          paste( levels( f )[ref[[s]]], collapse=collapse )
                          }
                        else names( ref )[s]
          newnames[ref[[s]]] <- rep( uninames[s], length( ref[[s]] ) )
        }
      levels( fnew ) <- newnames
      if( !is.null( first ) )
        {
      if( !first ) fnew <- factor( fnew, c( levels( f )[-unlist( ref )], uninames ) )
      if(  first ) fnew <- factor( fnew, c( uninames, levels( f )[-unlist( ref )] ) )
        }
    }

  # This is in order to merge levels with identical names
  #
  factor( fnew, levels=levels(fnew) )
}
