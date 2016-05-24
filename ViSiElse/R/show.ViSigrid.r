#' Method show for ViSigrid object.
#' @title Method \code{show-ViSigrid}
#' @name show-ViSigrid-method
#' @rdname show-ViSigrid-methods
#' @aliases show,ViSigrid-method
#' @exportMethod show
#' @docType methods
#' @param object a ViSigrid.
#' @seealso  \code{\linkS4class{ViSigrid}} and see \code{\link{plot-ViSigrid-method}} for examples.
setMethod( "show" , "ViSigrid" , function(object ) { 
  cat( "-parameters \n" )
  temp <- c()
  for (ii in seq( 1 , length( slot( object , "parameters" ) ) , 1 ) ) {
    temp <- paste(	temp, 
                   names(slot( object , "parameters" ))[ ii ] , 
                   paste( rep(" ",15 - nchar( names(slot( object , "parameters" ))[ ii ] ) ) ,   collapse =  "" ) ,
                   ": " ,
                   slot( object , "parameters" )[[ ii ]] , 
                   "\n"
    ) 
  }
  cat( temp )
  cat( "-MATp " ) 
  cat( paste( paste( rep(" ", 11 ) ,   collapse =  "" ), ":" ,
              dim( slot( object , "MATp" ) ) [ 1 ] ,"x", dim( slot( object , "MATp" ) ) [ 2 ] , "sparse Matrix of class \"dgCMatrix\" \n") )
  cat( "-L " ) 
  cat( paste( paste( rep(" ", 14 ) ,   collapse =  "" ),":" ,
              dim( slot( object , "L" ) ) [ 1 ] ,"x", dim( slot( object , "L" ) ) [ 2 ] , " data.frame \n") )
  cat( "-idsort " )
  cat( paste( paste( rep(" ", 9 ) ,   collapse =  "" ),":" ,
              dim( slot( object , "idsort" ) ) [ 1 ] ,"x", dim( slot( object , "idsort" ) ) [ 2 ]  , "matrix \n") )
  cat( "-MATpsup " ) 
  cat( paste( paste( rep(" ", 8 ) ,   collapse =  "" ),":" ,
              dim( slot( object , "MATpsup" ) ) [ 1 ] ,"x", dim( slot( object , "MATpsup" ) ) [ 2 ] , "sparse Matrix of class \"dgCMatrix\" \n") )
  cat( "-idsup " )
  cat( paste( paste( rep(" ", 10 ) ,   collapse =  "" ),":" ,
              "length " , length( slot( object , "idsup" ) )  , "vector \n") )
  cat( "-Lsup " ) 
  cat( paste( paste( rep(" ", 11 ) ,   collapse =  "" ),":" ,
              dim( slot( object , "Lsup" ) ) [ 1 ] ,"x", dim( slot( object , "Lsup" ) ) [ 2 ] , " data.frame \n") )
  cat( "-colvect " )
  cat( paste( paste( rep(" ", 8 ) ,   collapse =  "" ),":" ,
              dim( slot( object , "colvect" ) ) [ 1 ] ,"x", dim( slot( object , "colvect" ) ) [ 2 ]  , "matrix \n") )
  cat( "-BZL " ) 
  cat( paste( paste( rep(" ", 12 ) ,   collapse =  "" ),":" ,
              dim( slot( object , "BZL" ) ) [ 1 ] ,"x", dim( slot( object , "BZL" ) ) [ 2 ] , "sparse Matrix of class \"dgCMatrix\" \n") )
  cat( "-book " ) 
  cat( paste( paste( rep(" ", 11 ) ,   collapse =  "" ),":" ,
              dim( slot( object , "book" ) ) [ 1 ] ,"x", dim( slot( object , "book" ) ) [ 2 ] , " ViSibook \n") )
  cat( "-group " )
  cat( paste( paste( rep(" ", 10 ) ,   collapse =  "" ),":" ,
              "length " , length( slot( object , "group" ) ) , "factor \n") )
  cat( "-vect_tps " )
  cat( paste( paste( rep(" ", 7 ) ,   collapse =  "" ),":" ,
              "length " , length( slot( object , "vect_tps" ) ) , "vector \n") )
  cat( "-testsP " )
  cat( paste( paste( rep(" ", 9 ) ,   collapse =  "" ),":" ,
              "length " , length( slot( object , "testsP" ) )   , "vector \n") )
  cat( "-informers " )
  cat( paste( paste( rep(" ", 6 ) ,   collapse =  "" ),":" ,
              dim( slot( object , "informers" ) ) [ 1 ] ,"x", dim( slot( object , "informers" ) ) [ 2 ]  , "matrix \n") )
}
)
