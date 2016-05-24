detrend <-
function( M, t, weight=rep(1,nrow(M)) )
{
Thin.col <-
function ( X, tol = 1e-06)
# Function to remove lin. dep. columns from a matrix
{
  QR <- qr(X, tol = tol, LAPACK = FALSE)
  X[, QR$pivot[seq(length = QR$rank)], drop = FALSE]
}
# Now detrend the matrix using the weighted inner product.
  Thin.col( projection.ip( cbind( 1, t ), M , orth = TRUE, weight = weight ) )
}
