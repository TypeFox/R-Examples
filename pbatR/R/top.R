sortResults <- function( pbatObj, sortBy=NULL ) {
  ## I'm guessing Christoph titles this differently sometimes?
  if( is.null(sortBy) ) {
    guess <- c("powerFBAT")
  }else{
    guess <- sortBy
  }

  for( g in guess ) {
    wh <- which( g == names(pbatObj$results) )
    if( length(wh)>0 ) {
      pbatObj$results <- pbatObj$results[order(pbatObj$results[,wh], decreasing=TRUE),]
      return(pbatObj)
    }
  }

  stop( "pbat data could not be sorted. set the 'sortBy' option to the name that corresponds to the conditional power estimate you wish to use." )
}
top <- function( pbatObj, n=10, sortBy=NULL ) {
  if( class(pbatObj)[1]!="pbat" )
    stop( "Object must be of class 'pbat', i.e. a result of 'pbat.m(...)'." )

  pbatObj <- sortResults( pbatObj, sortBy=sortBy )

  if( n<1 || n>nrow(pbatObj$results) )
    n <- nrow( pbatObj$results )
  
  return( pbatObj$results[1:n,] )
}
