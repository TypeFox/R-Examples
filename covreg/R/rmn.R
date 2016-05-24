rmn <-
function (M = 0, Srow, Scol){
	m=dim(Srow)[1]
	n=dim(Scol)[1]
      tmp=eigen(Srow)
	Srow.h=tmp$vec %*% diag(sqrt(tmp$val),nrow=m) %*% t(tmp$vec)
	tmp=eigen(Scol)
	Scol.h=tmp$vec %*% diag(sqrt(tmp$val),nrow=n) %*% t(tmp$vec)
	Z=matrix(rnorm(m * n), m, n)
	Srow.h %*% Z %*% Scol.h + M
}
