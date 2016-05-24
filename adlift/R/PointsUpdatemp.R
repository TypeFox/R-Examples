"PointsUpdatemp" <-
function(X,coeff,nbrs,newnbrs,index,remove,pointsin,weights,lengths){

#does the update lifting step based on nbrs of remove

r<-which(pointsin==remove);

N<-length(pointsin)
pos<-NULL
for (i in 1:length(nbrs)){
	pos[i]<-min(which(newnbrs==nbrs[i]))
}

###update the interval lengths (>=2 nbrs)###

if ((r>=2)&(r<=(N-1))){
	lengths[index]<-as.row(lengths[index])
	weights<-as.row(weights)
	lengths[index]<-lengths[index]+lengths[r]*weights[pos]
}
else{
	if(r==1){
		lengths[2]<-lengths[2]+lengths[1]
	}
	if(r==N){
		lengths[N-1]<-lengths[N-1]+lengths[N]
	}
}  

###update the scaling function coefficients###

alpha<-matrix(0,1,length(nbrs))

if (length(nbrs)>=2){
	alpha<-lengths[r]*lengths[index]/(sum(lengths[index]^2))

	for (i in 1:length(nbrs)){
		coeff[[pointsin[index][i]]]<-coeff[[pointsin[index][i]]]+alpha[i]*coeff[[remove]]
	}
}
else{
	q<-which(pointsin==nbrs)
	alpha<-lengths[r]/lengths[q]
	coeff[[pointsin[q]]]<-coeff[[pointsin[q]]]+alpha*coeff[[remove]]
}

return(list(coeff=coeff,lengths=lengths,r=r,N=N,weights=weights,alpha=alpha))

}


