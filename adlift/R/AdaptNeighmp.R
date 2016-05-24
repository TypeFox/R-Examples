"AdaptNeighmp" <-function(pointsin,X,coefflist,coeff,nbrs,newnbrs,remove,intercept,neighbours,mpdet,g){

#does local adaptive prediction for the point remove based on N 
#points (chooses method of prediction and intercept);

mindetails<-NULL
minindices<-NULL
results<-list()
tempres<-list()
newinfo<-list()
nlist<-list()
clo<-NULL
int<-NULL
scheme<-NULL

N<-length(pointsin);
r<-which(pointsin==remove);

min1<-min(N-1,neighbours)
min2<-min(N-1,2*neighbours)

closest<-FALSE   #does prediction schemes with neighbours on each side

for (k in 1:min1){
	out1<-getnbrs(X,remove,pointsin,k,closest)
	nbrs<-out1$nbrs
	index<-out1$index
	newnbrs<-NULL
	for (i in 1:length(nbrs)){
		newnbrs<-c(newnbrs,rep(nbrs[i],times=length(g[[nbrs[i]]])))
	}
	out2<-AdaptPredmp(pointsin,X,coefflist,coeff,nbrs,newnbrs,remove,intercept,neighbours,mpdet,g)
	nlist[[k]]<-out1
	tempres[[k]]<-out2
}

closest<-TRUE

for (k in 1:min2){
	out1<-getnbrs(X,remove,pointsin,k,closest)
	nbrs<-out1$nbrs
	index<-out1$index
	newnbrs<-NULL
	for (i in 1:length(nbrs)){
		newnbrs<-c(newnbrs,rep(nbrs[i],times=length(g[[nbrs[i]]])))
	}
	out2<-AdaptPredmp(pointsin,X,coefflist,coeff,nbrs,newnbrs,remove,intercept,neighbours,mpdet,g)
	nlist[[k+min1]]<-out1
	tempres[[k+min1]]<-out2
}

for (i in 1:(min1+min2)){
	minindices[i]<-tempres[[i]][[7]]
	mindetails[i]<-tempres[[i]][[6]][minindices[i]]
}

totalminindex<-order(abs(mindetails))[1]

if(totalminindex<=neighbours){
	clo<-FALSE
}
else{
	clo<-TRUE
}

results<-tempres[[totalminindex]]
nbrs<-nlist[[totalminindex]]$nbrs
index<-nlist[[totalminindex]]$index

newinfo<-list(clo,nbrs,index)

return(list(results=results,newinfo=newinfo))
}


