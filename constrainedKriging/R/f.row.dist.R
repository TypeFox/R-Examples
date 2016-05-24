f.row.dist <-   function( matrix.1, matrix.2) {
    
    # purpose:          function computes the matrix of Euclidean distances between any pair 
    #                   of rows of the matrices matrix.1 and matrix.2. It returns the distance
    #                   matrix with dimensions c(nrow(matrix.1), nrow(matrix.2))
    # author:           A. Papritz
    # date:             2006-12-04
    
    if( !(is.matrix( matrix.1 ) & is.matrix( matrix.2 ) ) 
    )   stop( "the arguments must be matrices" )
    
    
    if( ncol( matrix.1 ) != ncol( matrix.2 ) 
    )   stop( "the matrices have inconsistent dimensions")
    
    if( nrow( matrix.2 ) <= nrow( matrix.1 ) ) {
        
        matrix.1 <- t( matrix.1 )
        t.dist.matrix <- matrix( 0, nrow = ncol( matrix.1 ), ncol = nrow( matrix.2 ) )
        
        for( t.i in 1:nrow( matrix.2 ) ) {
            
#             print( t.i )
            t.diff <- matrix.1 - matrix.2[t.i,]
            t.dist.matrix[, t.i] <- sqrt( colSums( t.diff * t.diff ) )
            
        }
        
    } else {
        
        matrix.2 <- t( matrix.2 )
        t.dist.matrix <- matrix( 0, nrow = ncol( matrix.2 ), ncol = nrow( matrix.1 ) )
        
        for( t.i in 1:nrow( matrix.1 ) ) {
            
#             print( t.i )
            t.diff <- matrix.2 - matrix.1[t.i,]
            t.dist.matrix[, t.i] <- sqrt( colSums( t.diff * t.diff ) )
            
        }
        
        t.dist.matrix <- t( t.dist.matrix )
    }
    
    
    
    return(t.dist.matrix)
    
}
