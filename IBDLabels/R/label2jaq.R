label2jaq <-
function( label, phased = TRUE){
   
    out <- allJaq()

    jaq <- numeric( length( label )) 
    sel <- ifelse( phased, 2,3 )
    for( i in 1:length(label)){
        jaq[i] <- out[ out[,1]== label[i], sel]
        }
        
    return( jaq ) 
}
