"UndoPointsUpdate" <-
function(X,coeff,nbrs,index,remove,r,N,pointsin,gamweights,lengths,lengthrem){

alpha<-matrix(0,1,length(nbrs))

if ((r>1)&(r<=N)) {
	alpha<-lengths[index]*lengthrem/(sum(lengths[index]^2))
	coeff[nbrs]<-coeff[nbrs]-alpha*coeff[remove]
	lengths[index]<-as.row(lengths[index])
	prod<-gamweights*lengthrem
	prod<-as.row(prod)
	lengths[index]<-lengths[index]-prod
}

if ((r==1)|(r==(N+1))) {
	q<-which(pointsin==nbrs)
	alpha<-lengthrem/lengths[q]
	coeff[pointsin[q]]<-coeff[pointsin[q]]-alpha*coeff[remove]
	lengths[q]<-lengths[q]-lengthrem
}

return(list(coeff=coeff,lengths=lengths,alpha=alpha))

}
