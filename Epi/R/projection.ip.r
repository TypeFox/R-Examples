projection.ip <-
function( X, M, orth = FALSE, weight=rep(1,nrow(X)) )
# Generate the projection of M on span(X) w.r.t the inner
# product <x|y>=sum( x*w*y).
# Avoids computing the entire projection matrix
#   X %*% inverse( X'WX ) %*% (XW)'  by first computing
#         inverse( X'WX ) %*% (XW)'M
# (which is (p x p) %*% (p x n) %*% (n x k), i.e. (p x k) )
# and then premultiplying X (n x p) hence avoiding making
# a n x n matrix underway (note that n is large, p is small).
# Note multiplication by W (diagional matrix) is done by
# vector multiplication using the recycling facility of R.
{
  if( nrow(X) != length(weight) )
      stop( "Dimension of space and length of weights differ!" )
  if( nrow(X) != nrow(M) )
      stop( "Dimension of space and rownumber of model matrix differ!" )
  Pp <- solve( crossprod( X * sqrt(weight) ), t( X * weight ) ) %*% M
  PM <- X %*% Pp
  if (orth) PM <- M - PM
  else PM
}
