HierarchWords <-
function (res,Table){
        CoR<-res
        if((is.data.frame(CoR))|(is.matrix(CoR))){
            Rcoor<-CoR
         }else{
	if(!is.null(res$row$coord))
		Rcoor<-res$row$coord
	if(!is.null(res$ind$coord))
		Rcoor<-res$ind$coord
       }
   DocTerm<-Table
CHC<-function(X, groups=NULL){
	
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
	if (is.null(groups)){
		groups<-list()
		for (i in 1:nrow(Sim)){
			groups[[i]]<-i
		}
	}
	distclust<-numeric()
	clust<-list()
	i<-1

	indice<-nrow(d)-1

	while(indice>0){
	
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
			
			clust[[i]]<-vector(mode="list",length=2)
			clust[[i]][[1]]<-groups[[minfc]]
			clust[[i]][[2]]<-groups[[maxfc]]

			rownames(Sim)[minfc]<-colnames(Sim)[minfc]<-rownames(d)[minfc]<-colnames(d)[minfc]<-rownames(Cont)[minfc]<-colnames(Cont)[minfc]<-paste(rownames(Sim)[minfc],"-",rownames(Sim)[maxfc])

			if (minfc!=1){
				Sim[minfc,minfc-1]<-Sim[minfc-1,minfc]<-0.5*Sim[minfc,minfc-1]+0.5*Sim[maxfc,minfc-1]-0.5*abs(Sim[minfc,minfc-1]-Sim[maxfc,minfc-1])
				d[minfc,minfc-1]<-d[minfc-1,minfc]<-0.5*d[minfc,minfc-1]+0.5*d[maxfc,minfc-1]+0.5*abs(d[minfc,minfc-1]-d[maxfc,minfc-1])
				Cont[minfc-1,minfc]<-Cont[minfc,minfc-1]<-1
			}
			if (maxfc!=nrow(SimCont)){
				Sim[maxfc+1,minfc]<-Sim[minfc,maxfc+1]<-0.5*Sim[minfc,maxfc+1]+0.5*Sim[maxfc,maxfc+1]-0.5*abs(Sim[minfc,maxfc+1]-Sim[maxfc,maxfc+1])
				d[maxfc+1,minfc]<-d[minfc,maxfc+1]<-0.5*d[minfc,maxfc+1]+0.5*d[maxfc,maxfc+1]+0.5*abs(d[minfc,maxfc+1]-d[maxfc,maxfc+1])
				Cont[maxfc+1,minfc]<-Cont[minfc,maxfc+1]<-1
			}
		
			groups[[minfc]]<-c(groups[[minfc]],groups[[maxfc]])
			groups<-groups[-maxfc]
			Sim<-Sim[-maxfc,-maxfc]
			d<-d[-maxfc,-maxfc]
			Cont<-Cont[-maxfc,-maxfc]
			i<-i+1
		
			SimCont<-Sim*Cont
			DistCont<-d*Cont
			indice<-indice-1
              }
			
	clust<-clust
           
	hc<-hclust(dist(X))
	hc$height<-distclust
	hc$order<-sort(hc$order)
	grups_blocs<-list()
	grups_blocs[[1]]<-rep(0,nrow(X))
	for (i in 1:(length(clust)-1)){
		grups_blocs[[i+1]]<-grups_blocs[[i]]
		grups_blocs[[i+1]][c(clust[[i]][[1]],clust[[i]][[2]])]<-i
	}

	for(i in 1:nrow(hc$merge)){
		if (length(clust[[i]][[1]])==1&length(clust[[i]][[2]])==1){
			hc$merge[i,1]<-(-clust[[i]][[1]])
			hc$merge[i,2]<-(-clust[[i]][[2]])
		}else{
			if (length(clust[[i]][[1]])==1){
				hc$merge[i,1]<-(-clust[[i]][[1]])
			}else{
				hc$merge[i,1]<-grups_blocs[[i]][clust[[i]][[1]][1]]
			}
		
			if (length(clust[[i]][[2]])==1){
				hc$merge[i,2]<-(-clust[[i]][[2]])
			}else{
				hc$merge[i,2]<-grups_blocs[[i]][clust[[i]][[2]][1]]
			}
		}
	}

	return (res=list(hc=hc, clust= clust, groups=groups, dist=distclust))
}

objCCC<-CHC(Rcoor)

espcron<-function(objCCC, DocTerm){
    mat<-DocTerm[which(rownames(DocTerm)%in%rownames(Rcoor)),]
	    FreqDoc<-apply(mat,2,sum)
	    mat<-mat[,which(FreqDoc>0)]
	
	nodos<-descfreq(mat)
      
	juntar<-function(elem){
		elem<-unlist(elem)
		aux<-as.data.frame(mat)
		aux[nrow(aux)+1,]<-apply(aux[elem,],2,sum)
		rownames(aux)[nrow(aux)]<-paste(rownames(aux)[min(elem)],"-",rownames(aux)[max(elem)],sep="")
		aux<-aux[-elem,]}
        numNodo<-vector("list",length(nodos))
	  for (i in 1:length(nodos)){
            Nodo<-as.vector(rep(i,length(nodos[[i]])))
            nodos[[i]]<- cbind(nodos[[i]],Nodo)
		}
         
                    
	pal.car<-lapply(lapply(objCCC$clust,function(el) juntar(el)),function(m) descfreq(m))
 	length(pal.car)
	length(pal.car[[1]])
       
     pal.carP<-pal.car
      numNodo<-vector("list",length(pal.carP))
       K <-length(pal.carP)+1
	for (i in 1:length(pal.carP)){
          for (j in 1:length(pal.carP[[i]])){
		Nodo<-as.vector(rep(i,length(pal.carP[[i]][[j]])))
           if(j<length(pal.carP[[i]])){   
           pal.carP[[i]][[j]]<-cbind(pal.carP[[i]][[j]],Nodo)
           }else{ 
             Nodo<-as.vector(rep(K+1,length(pal.carP[[i]][[j]])))
		pal.carP[[i]][[j]]<-cbind(pal.carP[[i]][[j]],Nodo)
		K=K+1 
          }
	   
 	   }  
        
       }
  
	res<-cbind(rep(names(nodos)[1],nrow(nodos[[1]])),nodos[[1]])
	for (i in 2:length(nodos)){
		res<-rbind(res,cbind(rep(names(nodos)[i],nrow(nodos[[i]])),nodos[[i]]))
	}

	for (i in 1:(length(pal.car)-1)){
		for (j in 1:length(pal.car[[i]])){
			res<-rbind(res,cbind(rep(names(pal.carP[[i]])[j],nrow(pal.carP[[i]][[j]])),pal.carP[[i]][[j]]))
		}	
	}

	res<-as.data.frame(res)
	colnames(res)[1]<-"Group"
	res$names<-rownames(res)
	res<-res[!duplicated(res),]
	for (i in 2:(ncol(res)-1)) res[,i]<-as.numeric(as.character(res[,i]))
	res[,1]<-as.character(res[,1])
	res<-res[which(res$v.test>0),]
	words<-levels(as.factor(res$names))
	char.words<-data.frame(matrix(nrow=length(words),ncol=ncol(res)-1))
	for (i in 1:length(words)){
		if (length(which(res$p.value==min(res$p.value[which(res$names==words[i]&res$v.test>0)])))==0) char.words[i,]<-0
		else char.words[i,]<-res[which(res$p.value==min(res$p.value[which(res$names==words[i]&res$v.test>0)])),1:8]
	}
	rownames(char.words)<-words
	colnames(char.words)<-colnames(res)[1:8]	
      char.words <- with(char.words, char.words[order(Nodo, decreasing = TRUE),])
      char.words<-char.words[,c(1,8,2:7)]
      colnames(char.words)[2] <- "Node"
	return(char.words)
 }	
  Hierarch.words<-espcron(objCCC, DocTerm)
  return(Hierarch.words)

}
