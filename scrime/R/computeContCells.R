`computeContCells` <-
function(data,computeExp=TRUE,justDiag=FALSE,check=TRUE,n.cat=NULL){
	if(!is.matrix(data))
		stop("data has to be a matrix.")
	if(is.null(n.cat))
		n.cat<-max(data,na.rm=TRUE)
	if(check)
		checkCatMat(data,n.cat)
	n.obs<-ncol(data)
	n.snp<-nrow(data)
	n.dist<-n.snp*(n.snp-1)/2
	listIdentity<-vector("list",n.cat)
	for(i in 1:n.cat)
		listIdentity[[i]]<-!is.na(data) & data==i
	mat.dist<-matrix(0,n.dist,ifelse(justDiag,n.cat,n.cat^2))
	n.add<-ifelse(justDiag,0,n.cat)
	for(i in 1:n.cat){
		tmp<-listIdentity[[i]]%*%t(listIdentity[[i]])
		mat.dist[,i+(i-1)*n.add]<-as.dist(tmp)
	}
	lowerIdentifier<-row(tmp)>col(tmp)
	if(!justDiag){
		upperIdentifier<-row(tmp)<col(tmp)
		id.up<-row2col(n.snp)
		for(i in 1:(n.cat-1)){
			for(j in (i+1):n.cat){
				tmp<-listIdentity[[i]]%*%t(listIdentity[[j]])
				mat.dist[id.up,(i-1)*n.cat+j]<-tmp[upperIdentifier]
				mat.dist[,(j-1)*n.cat+i]<-tmp[lowerIdentifier]
			}
		}
		tmp.names<-paste("N",rep(1:n.cat,e=n.cat),rep(1:n.cat,n.cat),sep="")
	}
	else
		tmp.names<-paste("N",1:n.cat,1:n.cat,sep="")
	colnames(mat.dist)<-tmp.names
	naIdentifier<-!is.na(data)
	n.obs<-naIdentifier%*%t(naIdentifier)
	if(computeExp){
		listRowSums<-lapply(listIdentity,function(x) x%*%t(naIdentifier))
			#lapply(listIdentity,rowSums,na.rm=TRUE)
		mat.exp<-matrix(0,n.dist,ifelse(justDiag,n.cat,n.cat^2))
		for(i in 1:n.cat){
			tmp<-(listRowSums[[i]]*t(listRowSums[[i]]))/n.obs
			mat.exp[,i+(i-1)*n.add]<-tmp[lowerIdentifier]
		}
		if(!justDiag){
			for(i in 1:(n.cat-1)){
				for(j in (i+1):n.cat){
					tmp<-(listRowSums[[i]]*t(listRowSums[[j]]))/n.obs
					mat.exp[id.up,(i-1)*n.cat+j]<-tmp[upperIdentifier]
					mat.exp[,(j-1)*n.cat+i]<-tmp[lowerIdentifier]
				}
			}
		}
		colnames(mat.exp)<-tmp.names
	}
	else
		mat.exp<-NULL
	return(list(mat.obs=mat.dist,mat.exp=mat.exp))
}

