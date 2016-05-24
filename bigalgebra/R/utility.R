is_transposed = function( tcode )
{
  if ( sum(tcode == c('n', 'N')) > 0 )
    return(FALSE)
  if ( sum(tcode == c('T', 't', 'C', 'c')) > 0 )
    return(TRUE)
  stop("Invalid transpose code given") 
}

# Throws an error if the matrix A is not among the listed classes
# and types, and returns TRUE if A is a big.matrix, FALSE otherwise.
check_matrix = function(A, classes=c('big.matrix', 'matrix'), 
  types='double')
{
  if (!any( class(A) == classes))
  {
    stop("A is not the correct class type")
  }
  if (!any(typeof(A) == types))
  {
    stop("The matrix type is not correct")
  }
  return( ifelse( class(A) == 'big.matrix', TRUE, FALSE ) )
}

# Create a big.matrix of the specified dimensions using the
# options(bigalgebra.temp_pattern) and options(bigalgebra.tempdir) naming
# convention. Such big.matrices are often used to store computed output.
# The resulting big.matrix has a finalizer that deletes the backing file
# when the R object is garbage collected.
# m: number of rows
# n: number of columns
# type: optional type, defaults to double
# val:  optional fill-in value of type 'type.'
anon_matrix = function(m, n, type, val=NULL)
{
  if(missing(type)) type = "double"
  f = basename(tempfile(pattern=options("bigalgebra.temp_pattern")[[1]]))
  p = options("bigalgebra.tempdir")[[1]]()
  d = sprintf("%s.desc",f)
  ans = filebacked.big.matrix(m, n, type, init=val, backingfile=f, backingpath=p,
                              descriptorfile=d)
  address = capture.output(print(ans@address))
  path = paste(p,f,sep="//")
  if(bigdebug()) warning(paste("Creating anonymous maitrx ",path))
  assign(address, path, envir=.bigalgebra_env)
  reg.finalizer(ans@address, finalize_anon_matrix, onexit=TRUE)
  ans
}

bigdebug = function()
{
  ans=tryCatch(options("bigalgebra.DEBUG")[[1]], error=function(e) FALSE)
  if(is.null(ans)) ans = FALSE
  ans
}

# Finalizer for anonymous big.matrices created by big algebra. This function
# deletes the corresponding backing file and descriptor file.
finalize_anon_matrix = function(e)
{
  address = capture.output(print(e))
  path    = get(address, envir=.bigalgebra_env)
  if(!is.null(path))
  {
    if(bigdebug()) warning(paste("Removing ",path))
    rm(list=address, envir=.bigalgebra_env)
    unlink(path)
    unlink(paste(path,"desc",sep="."))
  }
}
