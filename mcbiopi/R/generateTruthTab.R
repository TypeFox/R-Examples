generateTruthTab<-function(ltree){
	if(!is(ltree,"logregtree"))
		stop("ltree must be an object of class logregtree.")
	model.var<-ltree$trees[,3]
	model.var<-sort(model.var[model.var!=0])
	model.var<-model.var[!duplicated(model.var)]
	mat.perms<-getPerms(length(model.var))
	colnames(mat.perms)<-paste("X",model.var,sep="")
	mat.bin<-matrix(0,nrow(mat.perms),max(model.var))
	mat.bin[,model.var]<-mat.perms
	pred.out<-evalTree(ltree$trees, mat.bin)
	mat.truth<-cbind(mat.perms,outcome=pred.out)
	mat.truth
}

evalTree <- function(matTree, newBin){
	n.row <- nrow(matTree)
	mat <- matrix(1, ncol=n.row, nrow=nrow(newBin))
	for(i in 1:n.row){
		if(matTree[i,2] == 3){
			mat[,i] <- newBin[,matTree[i,3]]
			if(matTree[i,4] == 1)
				mat[,i] <- 1 - mat[,i]
		}
	}
	n.row <- floor(n.row/2)
	for(i in n.row:1){
		if(matTree[i,2]==1)
			mat[,i] <- mat[,2*i] * mat[,2*i+1]
		if(matTree[i,2]==2)
			mat[,i] <- mat[,2*i] + mat[,2*i+1]
	}
	mat[mat>1] <- 1
	if(sum(matTree != 0) == 0)
		mat[,1] <- 0
	mat[,1]
}


	