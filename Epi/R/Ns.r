#-------------------------------------------------------------------------------
# A wrapper for ns that automatically  takes the smallest and largest knots
# as boundary knots without further ado, but also allows cenering
# around a reference and detrending by means of projection
Ns <- function( x, ref = NULL, 
                    df = NULL,
                 knots = NULL,
             intercept = FALSE,
        Boundary.knots = NULL,
               detrend = FALSE ) 
{
## Check sensibility of arguments
if( !is.null(ref) )
  {
  if( !is.vector(ref) )
    stop( "Argument 'ref' must be a scalar, but it is a ", class(ref), "." )
  if( is.vector(ref) & length(ref)>1 )
    stop( "Argument 'ref' must be a scalar, but has length ", length(ref), "." )
  if( intercept )
    {  
    warning( "ref= specified, hence intercept=TRUE is ignored")
    intercept <- FALSE
    }
  }
## Detrending required?
if( any(detrend>0) ) ## covers both logical and vector
  {
  if( any(detrend<0) )
    stop( "Some elements of weight are <0, e.g. no",
          (ww <- which(detrend<0))[1:min(5,length(ww))], "." )
  if( !(length(detrend) %in% c(1,length(x))) )
    {
    warning( "Weights in inner product diagonal matrix set to 1")
    weight <- rep(1,length(x))    
    }
  else weight <- if( is.numeric(detrend) ) detrend else rep(1,length(x))
  detrend <- TRUE
  }
if( detrend & intercept )
  {  
  warning( "detrend= specified, hence intercept=TRUE is ignored")  
  intercept <- FALSE
  }
## Here is the specification of the spline basis
## df= specified
if( !is.null(df) )
  MM <- ns( x, df = df, intercept = (intercept & is.null(ref)) )
else
## knots= specified
{
if( is.null( Boundary.knots ) )
  {
  if( !is.null( knots ) )
    {
    knots <- sort( unique( knots ) )
    ok <- c(1,length(knots))
    Boundary.knots <- knots[ok]
    knots <- knots[-ok]
    }
  }
MM <- ns( x, knots = knots, Boundary.knots = Boundary.knots,
             intercept = (intercept & is.null(ref)) )
}
## Reference point specified ?
if( !is.null(ref) )
  {
  MM <- MM - ns( rep(ref,length(x)),
                 knots = attr(MM,"knots"),
        Boundary.knots = attr(MM,"Boundary.knots") )
  }
## Detrending required ?
if( detrend )
  { 
  DD <- detrend( MM, x, weight=weight )
  ## NOTE: detrend does not preserve attributes
  for( aa in c("degree","knots","Boundary.knots","intercept","class") )
     attr( DD, aa ) <- attr( MM, aa )    
  attr( DD, "detrend" ) <- TRUE
  attr( DD, "proj.wt" ) <- weight
  MM <- DD
  }
return( MM )
}
