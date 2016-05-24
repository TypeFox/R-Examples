#' Method print for ViSibook object.
#' @title Method \code{print-ViSibook}
#' @name print,ViSibook-method
#' @rdname print-ViSibook-methods
#' @aliases print,ViSibook-methods
#' print
#' print-ViSibook-methods
#' @exportMethod print
#' @docType methods
#' @param x a ViSibook object.
#' @seealso  \code{\linkS4class{ViSibook}},  \code{\link{buildViSiGrid}}, 
#' and see \code{\link{plot-ViSigrid-method}} for examples.
setMethod( "print", "ViSibook", function(x){ 
  tempcat <- paste0(	".......................\n" ,
                     " Book of the process : \n \n",
                     " >>Punctual actions : ", sum( methods::slot( x , "typeA" ) == rep( "p" , length( methods::slot( x , "vars" ) ) ) ) ,"\n"
  )
  # Punctuals Actions
  for (i in which( methods::slot( x , "typeA" ) == rep( "p" , length( methods::slot( x , "vars" ) ) ) ) ) {  # i=1
    # Name and label
    tempcat <- paste0( tempcat , "    -",methods::slot( x, "vars" ) [ i ] , " \n" , 
                       "        label : ",methods::slot( x, "label" ) [ i ] , " \n" ) 
    # Green Zone
    if ( is.na( match( "GZDeb" , methods::slot( x , "NAMES" ) ) ) == FALSE ) {
      if ( is.na( methods::slot( x, "GZDeb" ) [ i ]  ) == FALSE  ) {
        tempcat <- paste0( tempcat , 	"         Green zone   : " , methods::slot( x, "GZDeb" ) [ i ] , " to " , methods::slot( x, "GZFin" ) [ i ]," unit of time" ) 
        # Repetion 
        if ( is.na( methods::slot( x, "Repetition" ) [ i ]  ) == FALSE  ) {
          tempcat <- paste0( tempcat , ", repeated every ",methods::slot( x, "Repetition" ) [ i ] , " unit of time \n" ) 
        }else{tempcat <- paste0( tempcat , "\n") }
      }else{
        tempcat <- paste0( tempcat ,	"         Green zone   : NO  \n" ) 
      }
    }else{
      tempcat <- paste0( tempcat , "         Green zone   : NO  \n" ) 
    }
    # Black Zone 1
    if ( is.na( match( "BZBeforeDeb" , methods::slot( x , "NAMES" ) ) ) == FALSE ) {
      if ( is.na( methods::slot( x, "BZBeforeDeb" ) [ i ]  ) == FALSE  ) {
        tempcat <- paste0( tempcat , 	"         Black zone 1 : ", methods::slot( x, "BZBeforeDeb" ) [ i ] , " to " , methods::slot( x, "BZBeforeFin" ) [ i ] , " unit of time \n" ) 
      }else{
        tempcat <- paste0( tempcat , "         Black zone 1 : NO  \n" ) 
      }
    }else{
      tempcat <- paste0( tempcat , "         Black zone 1 : NO  \n" ) 
    }
    # Black Zone 2
    if ( is.na( match( "BZAfterDeb" , methods::slot( x , "NAMES" ) ) ) == FALSE ) {
      if ( is.na( methods::slot( x, "BZAfterDeb" ) [ i ]  ) == FALSE  ) {
        tempcat <- paste0( tempcat , 	"         Black zone 2 : ", methods::slot( x, "BZAfterDeb" ) [ i ] , " to " , methods::slot( x, "BZAfterFin" ) [ i ] , " unit of time \n " ) 
      }else{
        tempcat <- paste0( tempcat , "         Black zone 2 : NO  \n" ) 
      }
    }else{
      tempcat <- paste0( tempcat , 	"         Black zone 2 : NO  \n" ) 
    }
  }
  tempcat <- paste0( tempcat , 	"\n" , " >>Long actions : ", sum( methods::slot( x , "typeA" ) == rep( "l" , length( methods::slot( x , "vars" ) ) ) ) ," \n") 
  for (i in which( methods::slot( x , "typeA" ) == rep( "l" , length( methods::slot( x , "vars" ) ) ) ) ) {
    tempcat <- paste0( tempcat , 	"    -",methods::slot( x, "vars" ) [ i ] , " \n" , "       label : ",methods::slot( x , "label" ) [ i ] , " \n", "        Starts with : " , methods::slot( x , "deb" ) [ i ] , " \n",	"        Ends with   : " , methods::slot( x , "fin" ) [ i ] , " \n" ) 
    if ( is.na( match( "BZLong" , methods::slot( x , "NAMES" ) ) ) == FALSE ) {
      if ( is.na( methods::slot( x, "BZLong" ) [ i ]  ) == FALSE  ) {
        tempcat <- paste0( tempcat , 	"         Black zone : type ",slot( x, "BZLtype" ) [ i ]," -- ",slot( x, "BZLong" ) [ i ],  " unit of time \n" ) 
      }
    }else{
      tempcat <- paste0( tempcat , 	"         Black zone : NO  \n" ) 
    }
  }
  tempcat <- paste0( tempcat , "\n" , " >> Process to be drawn : ", sum( is.na( methods::slot( x , "showorder" ) ) == FALSE ) ,
                     " action", switch( as.character( sum( is.na( methods::slot( x , "showorder" ) ) == FALSE ) == 1) , "TRUE" = NA,"FALSE" = "s") , "\n") 
    sortindex <- sort( methods::slot( x , "showorder")  , index.return = TRUE)$ix
  if ( any( is.na( methods::slot( x , "showorder" ) ) ) ) {
    for (i in seq( 1 , sum( is.na( methods::slot( x , "showorder" ) ) ) , 1) ) { 
      sortindex <- mapply( FUN = function(x , y )(return( if ( y >= x ) { return( y + 1) }else{return( y ) } ) ) , x = which( is.na( methods::slot( x , "showorder" ) ) )[ i ] , y = sortindex )  
    }
  }
  for (i in sortindex ) {
    tempcat <- paste0( tempcat ,"          ", methods::slot( x , "showorder")[ i ] ," >> ",  methods::slot( x , "vars")[i]  )
    if (i != max( sortindex ) ) { tempcat <- paste0( tempcat ,"\n                   v \n" ) }
  }
  tempcat <- paste0( tempcat ,"\n \n .......................\n " )
  cat( tempcat ) 
}
)
