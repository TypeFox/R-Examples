"find.m" <-
function(
	M,	#matrix of a network
	clu,	#partition
	alt.blocks="reg", #alternative block to null block (for now only "reg" is supported)
	diag=!is.list(clu) ,#allow diagonal blocks
	cormet="none", #should we correct for diferent maxismum error contributins
			# "censor" - censor values larger than m
			# "correct" -  so that the maxsimum possible error contribution of the cell is the same regardles of a condition (either that somthing must be o or at least m)
	half = TRUE,	# should the returned value of m be one half of the value where the incosnistencies are the same, otherwise, the m is restricted to max(M)
	FUN="max"
){
	mx<-max(M)*(1+ half)
	mn<-min(M)
  diag=diag
	if(is.list(clu)){
		k<-sapply(clu,function(x)length(unique(x)))
		clu<-lapply(clu,function(x)as.integer(factor(x)))
		if(length(k)>2) {
			for(i in 2:length(clu)){
				clu[[i]]<-clu[[i]] + max(clu[[i-1]])
	  		}
	  		clu<-unlist(clu)
	  		clu<-list(clu,clu)
	  	}
	} else {
		clu<-as.integer(factor(clu))
		clu<-list(clu,clu)
		k<-sapply(clu,function(x)length(unique(x)))
	}
	
	m<-matrix(NA,nrow=k[1],ncol=k[2])
	err<-list(
		reg=function(B,m,FUN){
			nr<-dim(B)[1]	#numer of rows
			nc<-dim(B)[2]	#numer of colums
			sr<-apply(B,1,FUN);er<-m-sr[sr<m]
			sc<-apply(B,2,FUN);ec<-m-sc[sc<m]
			return(sum(er)*nc+ sum(ec)*nr - sum(pmin(rep(er,times=length(ec)),rep(ec,each=length(er)))))	#regular block error
		},
		com=function(B,m,FUN){sumpos(m-B)}
	)
	
	errd<-list(
		com=function(B,m,FUN){sumpos(m-B) + min(0,sum(diag(B))-sumpos(m-diag(B)))},
		reg=function(B,m,FUN){
			nr<-dim(B)[1]	#numer of rows
			nc<-dim(B)[2]	#numer of colums
			sr<-apply(B,1,FUN);er<-m-sr[sr<m]
			sc<-apply(B,2,FUN);ec<-m-sc[sc<m]
			return(sum(er)*nc+ sum(ec)*nr - sum(pmin(rep(er,times=length(ec)),rep(ec,each=length(er))))) #regular block error
		}
	)
	errd.null<-function(B,m){sum(B) - min(0,sum(diag(B))-sumpos(m-diag(B)))}
	
	for(i in 1:k[1]){
		for(j in 1:k[2]){
			B<-M[clu[[1]]==i,clu[[2]]==j, drop=FALSE]
			if(ss(B)==0) m[i,j]<-min(2*B[1,1],mx) else if(i==j&&diag&&sum(dim(B))>1){
				if(errd.null(B,m=mx)>=errd[[alt.blocks]](B,mx,FUN)*ifelse(cormet=="correct",(mx - 0)/(mx - mn),1)){
					m[i,j]<-mx
				}else{ 
					m[i,j]<-optimize(f=function(m,B,alt.blocks,FUN,cormet,mx,mn){corf<-ifelse(cormet=="correct", (mx - 0)/(m - mn),1); if(cormet=="censor") B[B>m]<-m;(errd.null(B,m)-errd[[alt.blocks]](B,m,FUN)*corf)^2},lower=ifelse(cormet=="censor",mn,0),upper=mx,B=B,FUN=FUN,alt.blocks=alt.blocks,cormet=cormet,mx=mx,mn=mn)$minimum
					if(cormet=="correct" && errd.null(B)<err[[alt.blocks]](B,m[i,j],FUN)*(mx - 0)/(m[i,j] - mn)) m[m<=min(B)]<-0				
			}
				
			}else{
				if(sum(B)>=err[[alt.blocks]](B,mx,FUN)*ifelse(cormet=="correct",(mx - 0)/(mx - mn),1)){
					m[i,j]<-mx 
				}else{
					m[i,j]<-optimize(f=function(m,B,alt.blocks,FUN,cormet,mx,mn){corf<-ifelse(cormet=="correct", (mx - 0)/(m - mn),1); if(cormet=="censor") B[B>m]<-m;(sum(B)-err[[alt.blocks]](B,m,FUN)*corf)^2},lower=ifelse(cormet=="censor",mn,0),upper=mx,B=B,FUN=FUN,alt.blocks=alt.blocks,cormet=cormet,mx=mx,mn=mn)$minimum
					if(cormet=="correct" && sum(B)<err[[alt.blocks]](B,m[i,j],FUN)*(mx - 0)/(m[i,j] - mn)) m[m<=min(B)]<-0				
				}
			}
		}
	}
	if(cormet=="censor") m[m<min(M[M>0])]<-0
	if(half) m<-m/2
	return(m)
}

