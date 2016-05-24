library( "linprog" )

# min x1 + x2, s.t. x1 + 0.5 * x2 = 2
cvec <- c( 1, 1 )
Amat <- matrix( c( 1, 0.5 ), nrow = 1 )
bvec <- 2
a1 <- solveLP( cvec, bvec, Amat, const.dir = "=" )
print( a1 )

a2 <- solveLP( cvec, bvec, Amat, const.dir = "=", lpSolve = TRUE )
print( a2 )

# max 27 * x1 + 9 * x2
# s.t. x1 - x2 = 8  &  x1 + x2 <= 74
cvec <- c( 27, 9 )
bvec <- c( 8, 74 )
Amat <- matrix( c( 1, 1, -1, 1 ), nrow = 2 ) 
b1 <- solveLP( cvec, bvec, Amat, maximum = TRUE, const.dir = c( "==", "<=" ) )
print( b1 )

b2 <- solveLP( cvec, bvec, Amat, maximum = TRUE, const.dir = c( "==", "<=" ),
   lpSolve = TRUE )
print( b2 )
