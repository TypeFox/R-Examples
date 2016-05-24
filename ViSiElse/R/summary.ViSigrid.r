#' Method summary for ViSigrid object.
#' @title Method \code{summary-ViSigrid}
#' @name summary-ViSigrid-method
#' @rdname summary-ViSigrid-methods
#' @aliases summary,ViSigrid-method
#' @exportMethod summary
#' @docType methods
#' @param object a ViSigrid.
#' @return list \itemize{
#' \item{ \strong{ punctuals} }{ summary of punctual actions (typeA=="p").}
#' \item{ \strong{ longs} }{  summary of long actions (typeA=="p"). }
#' }
#' @seealso \code{\linkS4class{ViSigrid}}, \code{\link{buildViSiGrid}},\code{\linkS4class{ViSibook}}.
#' and see \code{\link{plot-ViSigrid-method}} for examples.
setMethod( 	f = "summary" , 
            signature = "ViSigrid" , 
            definition = function(object ) (
              if ( is.null( methods::slot( object , "parameters")$informer  )) {
                cat( "No informers No tests were made in the call \n ")
              }else{ 
                cn = switch( methods::slot( object , "parameters")$informer , "median" = c( "q1","median","q3"), "mean" = c("mean-sd","mean","mean+sd" ) )
                infpunctuals <- methods::slot( object , "informers")[ , seq( 1 , sum( methods::slot( methods::slot(object , "book" ) , "typeA" ) == "p" ) , 1 ) ]
                rownames(infpunctuals) = rep(cn , dim(infpunctuals)[1]/3)
                if (length( methods::slot(object , "group" ) ) > 0 ) {
                  if (length( methods::slot( object , "testsP" ) ) > 0 ) { 
                    infpunctuals <- rbind( infpunctuals , methods::slot( object , "testsP")[ seq( 1 , sum( methods::slot( methods::slot(object , "book" ) , "typeA" ) == "p" ) , 1 ) ] )
                  }
                  rownames( infpunctuals ) <- c( paste(	rep("Gr" , 6) , 
                                                        c( rep( levels( methods::slot( object , "group" ) )[ 1 ] , 3 ) , 
                                                           rep( levels( methods::slot( object , "group" ) )[ 2 ] , 3 ) ) ,  
                                                        rep(cn , dim( infpunctuals )[ 1 ] / 3 ) ) , 
                                                 paste( switch( methods::slot( object , "parameters")$informer , "median" = "mood test", "mean" = "wilcoxon test" ) ," p.value < " , methods::slot( object , "parameters")$threshold.test) )
                }
                inflong <- methods::slot( object , "informers")[ , seq( sum( methods::slot( methods::slot(object , "book" ) , "typeA" ) == "p" ) + 1  , sum( methods::slot( methods::slot(object , "book" ) , "typeA" ) == "p" ) + sum( methods::slot( methods::slot(object , "book" ) , "typeA" ) == "l" ), 1 ) ]
                rownames(inflong) = rep(cn , dim(inflong)[1] / 3 )
                if (length( methods::slot(object , "group" ) ) > 0 ) {
                  if (length( methods::slot( object , "testsP" ) ) > 0 ) {
                    inflong <- rbind( inflong , methods::slot( object , "testsP")[ seq( sum( methods::slot( methods::slot(object , "book" ) , "typeA" ) == "p" ) + 1,  sum( methods::slot( methods::slot(object , "book" ) , "typeA" ) == "p" ) + sum( methods::slot( methods::slot(object , "book" ) , "typeA" ) == "l" ) , 1 ) ] )
                  }
                  rownames( inflong ) <- c( paste(	rep("Gr" , 6) , 
                                                   c( rep( levels( methods::slot( object , "group" ) )[ 1 ] , 3 ) , 
                                                      rep( levels( methods::slot( object , "group" ) )[ 2 ] , 3 ) ) ,  
                                                   rep(cn , dim( inflong )[ 1 ] / 3 ) ) , 
                                            paste( 	switch( methods::slot( object , "parameters")$informer , 
                                                            "median" = "mood test", 
                                                            "mean" = "wilcoxon test" ),
                                                    " p.value < " ,
                                                    methods::slot( object , "parameters")$threshold.test))				
                }
                return( list( punctuals = infpunctuals , longs = inflong ) )
              }
             )
)
