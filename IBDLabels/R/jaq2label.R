jaq2label <-
function( jaq, phased = TRUE){
    
    out <- allJaq()

    label <- numeric( length(jaq ))
    sel <- ifelse( phased, 2,3 )
    for( i in 1:length(jaq)) label[i] <- min( out[ out[,sel]==jaq[i], 1],
                                             na.rm = TRUE)


    if(!phased && length( intersect( jaq, c(3, 5, 7, 8 ))) > 0 )
        warning( "jaq9 goes to minimum label" ) 
    
    return( label )
}
