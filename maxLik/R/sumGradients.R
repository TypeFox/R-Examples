sumGradients <- function( gr, nParam ) {

   if( !is.null(dim(gr))) {
      gr <- colSums(gr)
   } else {
      ## ... or vector if only one parameter
      if( nParam == 1 && length( gr ) > 1 ) {
         gr <- sum(gr)
      }
   }
   return( gr )
}