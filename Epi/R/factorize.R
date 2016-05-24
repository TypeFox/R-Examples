# The factorize method
factorize <- function (x, ...) UseMethod("factorize")

# Default method is just the Relevel method
factorize.default <- Relevel.default

# The Lexis version of this
  Relevel.Lexis <-
factorize.Lexis <-
function (x, states=NULL, print=TRUE, ... )
{
# Is this really a Lexis object
if( !inherits(x,"Lexis") ) stop( "First argument must be a Lexis object" )

# If lex.Cst and lex.Xst are not factors, make them:
if( !is.factor(x$lex.Cst) | !is.factor(x$lex.Xst) )
  {
Cst <- factor(x$lex.Cst)
Xst <- factor(x$lex.Xst)
  }
else # just use the factors as they are
  {
Cst <- x$lex.Cst
Xst <- x$lex.Xst
  }
# - but amend them to have the same set of levels
all.levels = union(levels(Cst), levels(Xst))
Cst <- factor(Cst, levels = all.levels)
Xst <- factor(Xst, levels = all.levels)

# A table of actually occurring levels and their names
tCX <- table(Cst) + table(Xst)
all.levels <- names( tCX[tCX>0] )

# If states are not given, just return factors with reduced levels
if( is.null(states) )
  {
x$lex.Cst <- factor( Cst, levels = all.levels )
x$lex.Xst <- factor( Xst, levels = all.levels )
  }

# If new state names are given as a list it implies merging of them
if( !is.null( states ) & is.list( states ) )
  {
  x$lex.Cst <- Relevel( Cst, states, ... )
  x$lex.Xst <- Relevel( Xst, states, ... )
  if( print )
    {
    # Construct translation table between old and grouped states to print
    tC <- table( Cst, x$lex.Cst )
    tX <- table( Xst, x$lex.Xst )
    cC <- matrix( colnames(tC), nrow(tC), ncol(tC), byrow=T )
    cX <- matrix( colnames(tX), nrow(tX), ncol(tX), byrow=T )
    cC[tC==0] <- ""
    cX[tX==0] <- ""
    print( data.frame( type=rep( c("lex.Cst","lex.Xst"),
                                 c(nrow(tC),nrow(tX)) ),
                       old=c(rownames(tC),rownames(tX)),
                       new=c( apply( cC, 1, paste, collapse="" ),
                              apply( cX, 1, paste, collapse="" ) ) ) )
    }
  }

# If states is a character vector we assume that it's just new names
if( !is.null( states ) & is.character( states ) )
  {
  if( length( states ) != nlevels(Cst) )
    stop( "Second argument is a vector of length ", length(states),
          ", but it should be the joint no. of states, ",
         length(all.levels),
          "\ncorresponding to ", all.levels )
  levels( Cst ) <- levels( Xst ) <- states
  x$lex.Cst <- Cst
  x$lex.Xst <- Xst
  if( print )
    {
    cat( "New levels for lex.Xst and lex.Cst generated:\n" )
    print( data.frame( old=all.levels, new=levels(x$lex.Cst) ) )
    }
  }

# If states is a numeric vector we assume that it's just reordering
if( !is.null( states ) & is.numeric( states ) )
  {
  x$lex.Cst <- Relevel( Cst, states )
  x$lex.Xst <- Relevel( Xst, states )
  }

return( x )
}
