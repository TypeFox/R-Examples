
q0.xmin <- -10
q0.xmax <-   6
q0.n    <- 769
q0.dim  <-   5
q0.delx <- ( q0.xmax - q0.xmin )/( q0.n - 1 )

q0idx <- (function() {
    fudge <- ( ceiling( log( q0.n, base=2 ) ) 
               * .Machine$double.eps
               / q0.delx )
    
    function( x ) {
        floor( ( x - q0.xmin )/q0.delx + fudge ) + 1
    }
})()

q0sol <- function( x,i ) {
    x <- as.vector( x )
    y <- matrix( NA, length( i ), length( x ) )
    
    idx           <- q0idx( x )
    at.end        <- ( idx == q0.n )
    y[ ,at.end  ] <- q0.y[ i,q0.n ]

    x             <- x[ !at.end ]
    idx           <- idx[ !at.end ]
    y[ ,!at.end ] <- interp.h3( x, q0.x[ idx ],  q0.x[ idx+1 ],
                                   q0.y[ i,idx ],  q0.y[ i,idx+1 ],
                                   q0.yp[ i,idx ], q0.yp[ i,idx+1 ] )
    y
}
