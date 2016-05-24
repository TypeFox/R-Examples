uHomo.groups <-
function(res, N =1000,alfa=0.05){
	   CoR<-res
        if((is.data.frame(CoR))|(is.matrix(CoR))){
            Rcoor<-CoR
         }else{
	if(!is.null(res$row$coord))
		Rcoor<-res$row$coord
	if(!is.null(res$ind$coord))
		Rcoor<-res$ind$coord
       }

     X<- Rcoor

	### Similarity Matrix
	d<-dist(X)
	d0<-as.matrix(d)
	maxd<-max(d)
	maxd<-maxd+1e-10
	Sim<-as.matrix(maxd-d)
	Sim0<-Sim
	d<-as.matrix(d)
	
	### Constrained Matrix
	Cont<-matrix(nrow=nrow(Sim),ncol=ncol(Sim),0)
	Cont[1,2]<-1
	for (i in 2:(nrow(Cont)-1)){
		Cont[i,i+1]<-1
		Cont[i,i-1]<-1
	}
	Cont[nrow(Cont),nrow(Cont)-1]<-1
	rownames(Cont)<-rownames(Sim)
	colnames(Cont)<-colnames(Sim)

	### Similarity matrix used for constrained clustering
	SimCont<-Sim*Cont
	DistCont<-d*Cont
	
		groups<-list()
		for (i in 1:nrow(Sim)){
			groups[[i]]<-i
		}
	
	distclust<-numeric()
	clust<-list()
	i<-1

	### Clustering	
		while(sum(SimCont!=0)>0){
		
			### Find the position of the maxim similarity
			maxsim<-max(SimCont)
			posmaxsim<-which(SimCont==maxsim)
			
			if (posmaxsim[1]%%nrow(SimCont)==0){
				fila<-posmaxsim[1]%/%nrow(SimCont)
				col<-nrow(SimCont)
			}else{
				fila<-posmaxsim[1]%/%nrow(SimCont)+1
				col<-posmaxsim[1]%%nrow(SimCont)
			}
	
			maxfc<-max(fila,col)
			minfc<-min(fila,col)

			distclust[i]<-DistCont[fila,col]
		
			### Permutation test
			perm<-numeric()
			if (length(groups[[minfc]])==1&&length(groups[[maxfc]])==1){
				pvalid<-TRUE
			}else{
				##library(gdata)
				maux<-d0[c(groups[[minfc]],groups[[maxfc]]),c(groups[[minfc]],groups[[maxfc]])]	
				dmed<-median(upperTriangle(maux))
				mdist<-d0>dmed
				uns<-sum(mdist[groups[[minfc]],groups[[maxfc]]])
				ncomb<-factorial(length(groups[[minfc]])+length(groups[[maxfc]]))/(factorial(length(groups[[minfc]]))*factorial(length(groups[[maxfc]])))
				if (ncomb>N|ncomb==Inf|is.na(ncomb)){
					for (j in 1:N){
						elemperm<-sample(c(groups[[minfc]],groups[[maxfc]]),length(c(groups[[minfc]],groups[[maxfc]])))
						perm[j]<-sum(mdist[elemperm[1:length(groups[[minfc]])],elemperm[(length(groups[[minfc]])+1):length(elemperm)]])
					}
					if (sum(perm>=uns)<(N*alfa)) pvalid<-FALSE
					else pvalid<-TRUE
				}else{
					combin1<-combn(c(groups[[minfc]],groups[[maxfc]]),length(groups[[minfc]]),simplify=FALSE)
					combin2<-combn(c(groups[[minfc]],groups[[maxfc]]),length(groups[[maxfc]]),simplify=FALSE)

					for (j in 1:length(combin1)){
						perm[j]<-sum(mdist[combin1[[j]],combin2[[length(combin2)-j+1]]])
					}
	
					if (sum(perm>=uns)<(length(combin1)*alfa)) pvalid<-FALSE
					else pvalid<-TRUE
				}
			}

			### Join two clusters if permutation test is ok

			if (pvalid){		
				clust[[i]]<-vector(mode="list",length=2)
				clust[[i]][[1]]<-rownames(Sim)[minfc]
				clust[[i]][[2]]<-rownames(Sim)[maxfc]
				rownames(Sim)[minfc]<-colnames(Sim)[minfc]<-rownames(d)[minfc]<-colnames(d)[minfc]<-rownames(Cont)[minfc]<-colnames(Cont)[minfc]<-paste(rownames(Sim)[minfc],"-",rownames(Sim)[maxfc])

				if (minfc!=1){
					Sim[minfc,minfc-1]<-Sim[minfc-1,minfc]<-0.5*Sim[minfc,minfc-1]+0.5*Sim[maxfc,minfc-1]-0.5*abs(Sim[minfc,minfc-1]-Sim[maxfc,minfc-1])
					d[minfc,minfc-1]<-d[minfc-1,minfc]<-0.5*d[minfc,minfc-1]+0.5*d[maxfc,minfc-1]+0.5*abs(d[minfc,minfc-1]-d[maxfc,minfc-1])
				}
				if (maxfc!=nrow(SimCont)){
					Sim[maxfc+1,minfc]<-Sim[minfc,maxfc+1]<-0.5*Sim[minfc,maxfc+1]+0.5*Sim[maxfc,maxfc+1]-0.5*abs(Sim[minfc,maxfc+1]-Sim[maxfc,maxfc+1])
					d[maxfc+1,minfc]<-d[minfc,maxfc+1]<-0.5*d[minfc,maxfc+1]+0.5*d[maxfc,maxfc+1]+0.5*abs(d[minfc,maxfc+1]-d[maxfc,maxfc+1])
					Cont[maxfc+1,minfc]<-Cont[minfc,maxfc+1]<-Cont[maxfc,maxfc+1]
				}
				groups[[minfc]]<-c(groups[[minfc]],groups[[maxfc]])
				groups<-groups[-maxfc]
				Sim<-Sim[-maxfc,-maxfc]
				d<-d[-maxfc,-maxfc]
				Cont<-Cont[-maxfc,-maxfc]
				i<-i+1

			### If permutation test is not ok, change Constrained matrix 
		
			}else{
				Cont[maxfc,minfc]<-Cont[minfc,maxfc]<-0
			}

			SimCont<-Sim*Cont
			DistCont<-d*Cont
		} 
	return (res=list(NumHomogGroups=length(groups),groups=groups))
}
