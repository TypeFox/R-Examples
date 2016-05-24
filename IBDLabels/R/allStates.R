allStates <-
function( ngam ){

    nlabels <- maxlabel( ngam )
    out <- data.frame( label = 0:nlabels )

    out$vec <- allVec( ngam )
    out$lex <- allLex( ngam )

    if( ngam == 4 ){
        out$jaq15 <-  c(1,3,NA,4,2,5,6,9,11,10,7,13,12,14,8,15)
        out$jaq9 <- c(1,3,NA,3,2,4,5,7,8,7,5,8,8,8,6,9)
    }

    return( out) 
}
