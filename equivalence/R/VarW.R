"VarW" <-
function(W, piXXY, piXYY, n1, n2)
return( 1/(n1*n2) * (W - (n1+n2-1)*W^2 + (n1-1)*piXXY + (n2-1)*piXYY) )

