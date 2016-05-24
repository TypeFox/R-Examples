"invtnpmp" <-
function(X,coefflist,coeff,lengths,lengthsremove,pointsin,removelist,neighbrs,newneighbrs,schemehist,interhist,nadd=length(X)-2,intercept=TRUE,neighbours=1,closest=FALSE,LocalPred=LinearPredmp,mpdet="ave") {

xold<-NULL
fold<-NULL
g<-list()

xold<-X
fold<-coeff

isu<-adjustx(xold,xold,"mean")
X<-isu$sepx
g<-isu$groups

X<-as.row(X)
coeff<-as.row(coeff)
n<-length(X)
N<-length(pointsin)
m<-length(removelist)
d<-neighbours

nadd<-min(nadd,m)

for (i in 1:length(X)){
	coefflist[[i]]<-rep(coefflist[[i]],times=length(g[[i]]))
}


for (j in 1:nadd) {

	N<-length(pointsin)

	remove<-removelist[m-j+1]
	lengthrem<-lengthsremove[m-j+1]
	nbrs<-neighbrs[[m-j+1]]
	newnbrs<-newneighbrs[[m-j+1]]

	index<-NULL
	for (i in 1:length(nbrs)) {
		index[i]<-which(pointsin==nbrs[i])
	}

	B<-(X[remove]>X[nbrs])		#checks where to place removed point
	nt<-sum(B)			#by counting number of X positions on its 
					#left
	if (nt==0){
		r<-which(pointsin==nbrs[1])
	}
	if (nt==length(nbrs)){
		r<-which(pointsin==nbrs[length(nbrs)])+1
	}

	if ((nt>0)&(nt<length(nbrs))){			
		r<-which(pointsin==nbrs[nt+1])  
	} 

	if (is.null(schemehist)==FALSE){
		if (schemehist[m-j+1]=="Linear"){
       			res<-LinearPredmp(pointsin,X,coefflist,coeff,nbrs,newnbrs,remove,intercept=interhist[m-j+1],neighbours,mpdet,g)
		}
		if (schemehist[m-j+1]=="Quad"){
			res<-QuadPredmp(pointsin,X,coefflist,coeff,nbrs,newnbrs,remove,intercept=interhist[m-j+1],neighbours,mpdet,g)
		}
		if (schemehist[m-j+1]=="Cubic"){
			res<-CubicPredmp(pointsin,X,coefflist,coeff,nbrs,newnbrs,remove,intercept=interhist[m-j+1],neighbours,mpdet,g)
		}
	}
	else{ 
		res<-LocalPred(pointsin,X,coefflist,coeff,nbrs,newnbrs,remove,intercept,neighbours,mpdet,g)
	}
	if(length(res)==2){
		l<-res[[1]]
	}
	else{
		l<-res
	}
	gamweights<-l[[1]]

	l1<-UndoPointsUpdatemp(X,coeff,nbrs,newnbrs,index,remove,r,N,pointsin,gamweights,lengths,lengthrem)
	coeff<-l1$coeff
	lengths<-l1$lengths

	for (i in 1:length(nbrs)){
		coefflist[[nbrs[i]]]<-rep(coeff[nbrs[i]],length=length(g[[nbrs[i]]]))
	}

	if(length(nbrs)>1){
		coeff1<-NULL
		for (i in 1:length(nbrs)){
			coeff1<-cbind(coeff1,as.row(coefflist[[nbrs[i]]]))
		}
		pred<-sum(as.column(gamweights)*as.column(coeff1))
	}
	else{
		pred<-coeff[nbrs]
	}

	coeff[remove]<-coeff[remove]+pred
	coefflist[[remove]]<-rep(coeff[remove],length=length(g[[remove]]))

	removelist<-setdiff(removelist,remove)

	if (r==1){
		lengths<-c(lengthrem,lengths)
		pointsin<-c(remove,pointsin)
	}
	if (r==(N+1)){
		lengths<-c(lengths,lengthrem)
		pointsin<-c(pointsin,remove)
	}
	if ((r>1)&(r<(N+1))){
		lengths<-c(lengths[1:(r-1)],lengthrem,lengths[r:N])
		pointsin<-c(pointsin[1:(r-1)],remove,pointsin[r:N])
	}

}

return(list(X=X,coeff=coeff,lengths=lengths,lengthsremove=lengthsremove,pointsin=pointsin,removelist=removelist))
}

