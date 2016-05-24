Bsublin <-
function(X,R,R0=0)
{
	Bx = X
	Bx[X>R] = (R+R0)/2
	ind = X>R0 & X<=R
	Bx[ind] =(-X[ind]^2 +2*R*X[ind]-R0^2)/(2*(R-R0))
	return(Bx)
}
