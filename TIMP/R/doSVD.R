"doSVD" <-
function (x, numleft, numright)
#doSVD
{
	svdx<-La.svd(as.matrix(x),nu=numleft,nv=numright)
	list(values=svdx$d, left=as.matrix(svdx$u),
             right=as.matrix(svdx$vt))

}

