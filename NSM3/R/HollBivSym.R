HollBivSym<-function(x,y=NA){
	check<-0
	if((is.null(ncol(x))||ncol(x)==1)&&!is.na(y)){
		check=1
	}

	if(max(dim(x)[2],1)==2){
		y<-x[!is.na(x[,2]),2]
		x<-x[!is.na(x[,1]),1]
		check=1
	}
 
	if(!check){
		return('Error: invalid form for entered data')
	}

	obs.data<-cbind(x,y)
	a.vec<-apply(obs.data,1,min)
	b.vec<-apply(obs.data,1,max)

	n<-length(a.vec)
	test<-function(r,c) {as.numeric((a.vec[c]<b.vec[r])&&(b.vec[r]<=b.vec[c])&&(a.vec[r]<=a.vec[c]))}
	myVecFun <- Vectorize(test,vectorize.args = c('r','c')) 
	d.mat<-outer(1:n, 1:n, FUN=myVecFun) 

	A.calc<-function(r.vec){
		s.vec<-2*r.vec-1
		T.vec<-s.vec%*%d.mat
		A.obs<-sum(T.vec*T.vec)/n^2
		return(A.obs)
	}

	A.obs<-A.calc(apply(obs.data,1,function(x){x[1]<x[2]}))
	return(A.obs)
}



