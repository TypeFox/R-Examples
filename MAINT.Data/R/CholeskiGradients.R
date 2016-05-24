fL2 <- function(fx,L2,...)
{
	n = nrow(L2)
	p = ncol(L2)
	res = 0.
	for (i in 1:n) for (j in 1:p) res = res + fx(j1=i,j2=j,L2=L2,...)
	return(res)
}

fL <- function(fx,L=NULL,L2=NULL,...) 
{
	if (is.null(L) && is.null(L2)) stop("At leat one of the L, L2 arguments must be supplied\n")
	if (is.null(L2)) return(fL2(fx=fx,L2=L%*%t(L),...)) 
	if (!is.null(L2)) return(fL2(fx=fx,L2=L2,...))
}

fL.grad <- function(dfx,L,l1,l2,totald=FALSE,L2=NULL,...)
{
	p = ncol(L)
	colgrad = array(dim=p-l2+1)
	if (is.null(L2)) L2 = L%*%t(L)
	if (l1 < l2) { tmp=l1; l1=l2; l2=tmp }   # Force arguments to be consistent with a Low-Triangular matrix structure   
	for (a in l2:p)  {
		if (totald==FALSE) colgrad[a-l1+1] = dfx(j1=a,j2=l2,L2=L2,...) + dfx(j1=l2,j2=a,L2=L2,...)
		if (totald==TRUE) 
			if (a==l1) colgrad[a-l2+1] = 2*dfx(j1=l1,j2=l1,L2=L2,...)
			else colgrad[a-l2+1] = dfx(j1=max(a,l1),j2=min(a,l1),L2=L2,...)
	}
	return(colgrad%*%L[l2:p,l2])
}

