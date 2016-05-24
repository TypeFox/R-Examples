
########################################
########################################
####Functions for pre-processing data, interactions, etc.
########################################
########################################
make.inter<-function(X,treat){

	if(length(colnames(X))==0) colnames(X)<-paste("X",1:ncol(X),sep="")
	
	X<-X[,apply(X,2,sd)>0]
	n<-nrow(X)

	if(length(colnames(X))==0) colnames(X)<-paste("X_",1:ncol(X),sep="")

	treat.mat<-matrix(NA,nrow=n,ncol=sum(apply(treat, 2, FUN=function(x) length(unique(x)))))
	colnames(treat.mat)<-paste("dum",1:ncol(treat.mat))
	treat<-data.frame(treat)
	put.col<-1

	for(i.treat in 1:ncol(treat)){
		for(j.treat in sort(unique(treat[,i.treat]))){
			treat.mat[,put.col]<-1*(treat[,i.treat]==j.treat)
			colnames(treat.mat)[put.col]<-paste(colnames(treat)[i.treat],j.treat,sep="_")
			put.col<-put.col+1
		}
	}

	put.vec<-1
	names.vec<-NULL
	clust.vec<-NULL
	for(i.treat in 1:ncol(treat.mat)){
		for(i.X in 1:ncol(X)){
		names.vec[put.vec]<-paste(colnames(X)[i.X], colnames(treat.mat)[i.treat], sep=" x ")
		clust.vec[put.vec]<-i.treat
		put.vec<-put.vec+1
		}
	}
	
	names.vec<-gsub(" x int","",names.vec)
	names.vec<-gsub("int x ","",names.vec)

	X.temp<-matrix(NA,nrow=n, ncol=ncol(treat.mat)*ncol(X))

	big.X<- cbind(X[,apply(X,2,sd)>0], treat.mat, makeinter_cpp("X00"=treat.mat,"X10"=X,"X20" = X.temp)$X)
	colnames(big.X)<-c(colnames(X[,apply(X,2,sd)>0]),colnames(treat.mat),names.vec)
	clust.vec<-c(rep(1,ncol(X)),rep(2,ncol(treat.mat)),clust.vec+2)
return(list("X.big"=big.X,"c.clust"=clust.vec))
}


make.inter.treat<-function(X,treat){
	
	if(length(colnames(X))==0) colnames(X)<-paste("X",1:ncol(X),sep="")
	treat<-as.matrix(treat)
	
	n<-nrow(treat)

	treat.mat<-matrix(NA,nrow=n,ncol=sum(apply(treat, 2, FUN=function(x) length(unique(x)))))
	colnames(treat.mat)<-paste("dum",1:ncol(treat.mat))
	treat<-data.frame(treat)
	put.col<-1

	for(i.treat in 1:ncol(treat)){
		for(j.treat in sort(unique(treat[,i.treat]))){
			treat.mat[,put.col]<-1*(treat[,i.treat]==j.treat)
			colnames(treat.mat)[put.col]<-paste(colnames(treat)[i.treat],j.treat,sep="_")
			put.col<-put.col+1
		}
	}
	

	put.vec<-1
	names.vec<-NULL
	clust.vec<-NULL
	for(i.treat in 1:ncol(treat.mat)){
		for(i.treat2 in 1:ncol(treat.mat)){
		names.vec[put.vec]<-paste(colnames(treat.mat)[i.treat2], colnames(treat.mat)[i.treat], sep=" x ")
		clust.vec[put.vec]<-i.treat
		put.vec<-put.vec+1
		}
	}
	
	names.vec<-gsub(" x int","",names.vec)
	names.vec<-gsub("int x ","",names.vec)

	X.temp<-matrix(NA,nrow=n, ncol=ncol(treat.mat)^2)

	big.X<- cbind(X[,apply(X,2,sd)>0], treat.mat, makeinter_cpp("X00"=treat.mat,"X10"=treat.mat,"X20" = X.temp)$X)
	colnames(big.X)<-c(colnames(X[,apply(X,2,sd)>0]),colnames(treat.mat),names.vec)
	clust.vec<-c(rep(1,ncol(X)),rep(2,ncol(treat.mat)),clust.vec+2)
return(list("X.big"=big.X,"c.clust"=clust.vec))
}





make.threewayinter<-function(X,treat,treat2){
	n<-nrow(X)
	if(length(colnames(X))==0) colnames(X)<-paste("X",1:ncol(X),sep="")
	X<-X[,apply(X,2,sd)>0]
	X<-cbind(1,X)
	colnames(X)[1]<-"int"

	treat.mat<-matrix(NA,nrow=n,ncol=sum(apply(as.matrix(treat), 2, FUN=function(x) length(unique(x)))))
	colnames(treat.mat)<-paste("dum",1:ncol(treat.mat))
	treat<-data.frame(treat)
	put.col<-1

	for(i.treat in 1:ncol(treat)){
		for(j.treat in sort(unique(treat[,i.treat]))){
			treat.mat[,put.col]<-1*(treat[,i.treat]==j.treat)
			colnames(treat.mat)[put.col]<-paste(colnames(treat)[i.treat],j.treat,sep="_")
			put.col<-put.col+1
		}
	}

	treat.mat2<-matrix(NA,nrow=n,ncol=sum(apply(as.matrix(treat2), 2, FUN=function(x) length(unique(x)))))
	colnames(treat.mat2)<-paste("dum",1:ncol(treat.mat2))
	treat2<-data.frame(treat2)
	put.col<-1

	for(i.treat in 1:ncol(treat2)){
		for(j.treat in sort(unique(treat2[,i.treat]))){
			treat.mat2[,put.col]<-1*(treat2[,i.treat]==j.treat)
			colnames(treat.mat2)[put.col]<-paste(colnames(treat2)[i.treat],j.treat,sep="_")
			put.col<-put.col+1
		}
	}
	
	

	put.vec<-1
	names.vec<-NULL
	clust.vec<-NULL
	for(i.treat2 in 1:ncol(treat.mat2)){
	for(i.treat in 1:ncol(treat.mat)){
		for(i.X in 1:ncol(X)){
		names.vec[put.vec]<-paste(colnames(X)[i.X], colnames(treat.mat)[i.treat],colnames(treat.mat2)[i.treat2], sep=" x ")
		clust.vec[put.vec]<-i.treat*100000+i.treat2
		put.vec<-put.vec+1
		}
	}
	}	
	names.vec<-gsub(" x int","",names.vec)
	names.vec<-gsub("int x ","",names.vec)

	X.temp<-matrix(NA,nrow=n,ncol=length(names.vec))
	X.eval<-makethreeinter_cpp("X00"=treat.mat,"X10"=treat.mat,"X110"=X,"X20" = X.temp)
	big.X2<- X.eval$X
	dropcols<-is.na(big.X2[1,])|(colSums(big.X2^2,na.rm=TRUE)<.0000000000001)
	big.X2<-big.X2[,!dropcols]
	clustvec2<-clust.vec[!dropcols]
	colnames(big.X2)<-names.vec[X.eval$colvec[X.eval$colvec!=0]]
	rm(X.eval)
	
	##Make twoway interactions
	X0.eval<-make.inter(treat.mat,treat.mat)
	big.X0<-X0.eval$X
	dropcols<-colSums(big.X0^2)<.0000000000001
	dropcols[grep("_0",colnames(big.X0))]<-TRUE
	big.X0<-big.X0[,!dropcols]
	clustvec0<-X0.eval$c.clust[!dropcols]
	rm(X0.eval)
	
	X1.eval<-make.inter(X,treat)
	big.X1<-X1.eval$X
	clustvec1<-X1.eval$c.clust+max(clustvec0)
	rm(X1.eval)
	
		clustvec2[clustvec2>-1]<-max(c(clustvec1,clustvec2))+1
	clust.vec.try<-c(clustvec0,clustvec1,clustvec2)
	clust.vec.try<-as.numeric(as.factor(clust.vec.try))
return(list("big.X"=cbind(big.X0,big.X1,big.X2),"c.clust"=clust.vec.try))
}



make.justthreeway<-function(X1,X2,X3){
	n<-nrow(X1)
	add.int<-function(x){
		if(length(colnames(x))==0) colnames(x)<-paste("X",1:ncol(x),sep="")
		x<-x[,apply(x,2,sd)>0]
		#x<-cbind(x)
		#colnames(x)[1]<-"int"
		#x
	}
	
	X1<-add.int(X1)
	X2<-add.int(X2)
	X3<-add.int(X3)

	put.vec<-1
	names.vec<-NULL
	clust.vec<-NULL
	for(i.X3 in 1:ncol(X3)){
	for(i.X2 in 1:ncol(X2)){
		for(i.X1 in 1:ncol(X1)){
		names.vec[put.vec]<-paste(colnames(X1)[i.X1], colnames(X2)[i.X2],colnames(X3)[i.X3], sep=" x ")
		put.vec<-put.vec+1
		}
	}
	}	
	names.vec<-gsub(" x int","",names.vec)
	names.vec<-gsub("int x ","",names.vec)

	X.temp<-matrix(NA,nrow=n,ncol=length(names.vec))
	X.eval<-makethreeinter_cpp("X00"=X3,"X10"=X2,"X110"=X1,"X20" = X.temp)
	big.X2<- X.eval$X
	colnames(big.X2)<-names.vec
	dropcols<-is.na(big.X2[1,])|(colSums(big.X2^2,na.rm=TRUE)<.0000000000001)
	big.X2<-big.X2[,!dropcols]
	big.X2
}##Close out make.justthreeway


make.interXX<-function(X,X2){

	if(length(colnames(X))==0) colnames(X)<-paste("X",1:ncol(X),sep="")
	if(length(colnames(X2))==0) colnames(X2)<-paste("X2",1:ncol(X),sep="")
	
	X<-X[,apply(X,2,sd)>0]
	n<-nrow(X)


	treat.mat<-X2
	put.vec<-1
	names.vec<-NULL
	clust.vec<-NULL
	for(i.treat in 1:ncol(treat.mat)){
		for(i.X in 1:ncol(X)){
		names.vec[put.vec]<-paste(colnames(X)[i.X], colnames(treat.mat)[i.treat], sep=" x ")
		clust.vec[put.vec]<-i.treat
		put.vec<-put.vec+1
		}
	}
	
	names.vec<-gsub(" x int","",names.vec)
	names.vec<-gsub("int x ","",names.vec)

	X.temp<-matrix(NA,nrow=n, ncol=ncol(treat.mat)*ncol(X))

	big.X<- cbind(X[,apply(X,2,sd)>0], treat.mat, makeinter_cpp("X00"=treat.mat,"X10"=X,"X20" = X.temp)$X)
	colnames(big.X)<-c(colnames(X[,apply(X,2,sd)>0]),colnames(treat.mat),names.vec)
	clust.vec<-c(rep(1,ncol(X)),rep(2,ncol(treat.mat)),clust.vec+2)
return(list("X.big"=big.X,"c.clust"=clust.vec))
}





