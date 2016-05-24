allLex <-function( ngam ){

    nlabels <- maxlabel( ngam )
    invalid <- sapply( 0:nlabels, function(x){
                          any(is.na(label2vec( x, ngam = ngam)))})
 
    all.lex <- rep( NA, length = nlabels+1 )
    all.lex[ !invalid ] <- 1:sum( !invalid )
    names( all.lex ) <- 0:nlabels
    return( all.lex )
}
    
    
######### OLD VERSION ##########################################

## This one made up all the vectors and stored them - just need to
## know if it was a warning or not

## allLex <-
## function( ngam ){
##     ## All lexicographic labels

##     ## library(mgcv) ## for uniquecombs

##     nlabels <- maxlabel( ngam )
##     all.vecs <- allVec( ngam )
##     all.simp <- t( apply( all.vecs, 1, fgl2label )  )
##     all.lex <- attr( uniquecombs( all.simp ), "index" )

##     names( all.lex ) <- 0:nlabels
##     return( all.lex )
## }
