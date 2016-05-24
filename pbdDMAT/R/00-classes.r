# --------------------------------------------------
# Validity methods
# --------------------------------------------------

# Validity checking for ddmatrix objects
valid.ddmatrix <- function(object)
{
  # check for valid context
  test <- exists(paste(".__blacs_gridinfo_", object@ICTXT, sep=""), envir=.pbdBASEEnv)
  if (!test)
    return(paste("Context", object@ICTXT, "is not a valid context"))
  
  
  # check that the dims are theoretically reasonable
  if ( !(is.numeric(object@dim) && length(object@dim)==2) )
    return("Invalid slot 'dim'")
  if ( !(is.numeric(object@ldim) && length(object@ldim)==2) )
    return("Invalid slot 'ldim'")
  if ( !(is.numeric(object@bldim) && length(object@bldim)==2) )
    return("Invalid slot 'bldim'")
  
  
  # check valid ldim (assuming valid dim, ldim, and ictxt)
  ldim <- base.numroc(dim=object@dim, bldim=object@bldim, ICTXT=object@ICTXT, fixme=TRUE)
  
  if ( !all(ldim==dim(object@Data)) )
    return("dim(Data) not valid for this choice of 'dim', 'bldim', and 'ICTXT'")
  
  # undable to find a problem...
  return(TRUE)
}

# --------------------------------------------------
# Matrices
# --------------------------------------------------

setClassUnion("Linalg", c("vector", "matrix"))

setClass(Class="dmat", 
  representation=representation(
    Data="Linalg", 
    dim="numeric", 
    ldim="numeric",
    storage="character"#, "VIRTUAL"
    ),
  prototype=prototype(
#    Data=matrix(0.0), 
    dim=c(1L, 1L), 
    ldim=c(1L, 1L),
    storage="llb") # locally load balanced
)



#' Class ddmatrix
#' 
#' Distributed matrix class.
#' 
#' @slot DATA
#' The local submatrix.
#' @slot bldim
#' Blocking factor.
#' @slot ICTXT
#' BLACS ICTXT value.  Should be one of 0, 1, or 2 (initialized from
#' \code{pbdBASE::init.grid()}) or a custom value greater than 2 (created from
#' \code{pbdBASE::blacs_gridinit()}).
#' 
#' @name ddmatrix-class
#' @keywords Classes
#' @docType class
setClass(
  Class="ddmatrix", 
  representation=representation(
    Data="matrix",
    bldim="numeric",
    ICTXT="numeric"
  ),
  
  prototype=prototype(
    Data=matrix(0.0),
    dim=c(1L, 1L),
    ldim=c(1L, 1L),
    bldim=c(1L, 1L),
    ICTXT=0L,
    storage="scalapack"
  ),
  contains="dmat"#,
  #
  #validity=valid.ddmatrix
)

# Distributed Sparse Matrix
setClass(
  Class="dsmatrix", 
  representation=representation(
    Data="numeric",
    row_ptr="numeric",
    col_ind="numeric"
    
  ),
  
  prototype=prototype(
    Data=0.0,
    dim=c(1L, 1L),
    ldim=c(1L, 1L),
    row_ptr=1,
    col_ind=1,
    storage="csr"
  ),
  contains="dmat",
)



# --------------------------------------------------
# Vectors
# --------------------------------------------------

# Distributed Dense Vector
setClass(
  Class="ddvector", 
  representation=representation(
     Data="vector",
     len="numeric",
     llen="numeric",
     bldim="numeric",
     ICTXT="numeric"
  ),
  
  prototype=prototype(
    Data=0.0,
    len=1L,
    llen=1L,
    bldim=c(1L, 1L),
    ICTXT=0L
  )
)


# Distributed Sparse Vector
setClass(
  Class="dsvector", 
  representation=representation(
    Data="vector",
    length="numeric",
    llength="numeric",
    row_ptr="numeric",
    col_ind="numeric",
    storage="character"
  ),
  
  prototype=prototype(
    Data=0.0,
    length=1L,
    llength=1L,
    row_ptr=1,
    col_ind=1,
    storage="csr"
  )
)
