`replaceNAs` <-
function(mat,mat.na,nn=3,distance="smc2Mats",n.num=100,use.weights=TRUE,n.cat=NULL){
	n.row<-nrow(mat.na)
	FUN<-match.fun(distance)
	sets<-c(seq(1,n.row,n.num),n.row+1)
	sets<-unique(sets)
	for(i in 1:(length(sets)-1)){
		consider<-sets[i]:(sets[i+1]-1)
		tmp<-mat.na[consider,,drop=FALSE]
		mat.dist<-FUN(mat,tmp,n.cat=n.cat)
		colS<-colSums(mat.dist < 10^-8)
		ids1<-which(colS>0)
		ids2<-which(colS==0)
		for(j in ids1){
			tmp.ids<-which(mat.dist[,j] < 10^-8)[1]
			tmp[j,]<-mat[tmp.ids,]
		}
		for(j in ids2){
			tmp.dist<-mat.dist[,j]
			names(tmp.dist)<-1:length(tmp.dist)
			tmp.dist<-sort(tmp.dist)[1:nn]
			tmp.ids<-as.numeric(names(tmp.dist))
			tmp.mat<-mat[tmp.ids,,drop=FALSE]
			tmp.ids2<-which(is.na(tmp[j,]))
			for(k in tmp.ids2)
				tmp[j,k]<-if(use.weights) weightMode(tmp.dist^-1,tmp.mat[,k])
					else modeDist(tmp.mat[,k])
		}
		mat.na[consider,]<-tmp
	}
	mat.na
}

