ordi.breadth <-
function(dat,dist.method="jaccard"){#dat is data with herbivore species as rows and plants as column
	
				dat.by.host<-t(dat)
				dismatrix<-vegdist(dat.by.host,method=dist.method,binary=TRUE)
				distances<-vector("list",dim(dat)[1])
			#index.for.breadth <- 1:dim(dat)[1]
			host.breadth<-NA
			tot.breadth<-numeric(length=dim(dat)[1])
			centroids.group<-array(dim=c(dim(dat)[1],dim(dat)[2]-1))
						#define group vector
						group.vectors<-array("NO",dim=c(dim(dat)))
						group.vectors[which(dat!=0)]<-"YES"
		
		ug<-betadisperF(dismatrix,group=rep("YES",dim(dat)[2]),type="centroid")
		
		index.for.breadth<-1:dim(dat)[1]
		for(i in index.for.breadth){
		tempdisper<-betadisperF(dismatrix,group=group.vectors[i,],type="centroid")
		distances[[i]]<-tempdisper$distances
		tot.breadth[i]<-sum(tempdisper$distances[which(group.vectors[i,]=="YES")])
		centroids.group[i,1:length(tempdisper$centroids[rownames(tempdisper$centroids)=="YES",])]<-tempdisper$centroids[rownames(tempdisper$centroids)=="YES",]			
					}
		scaled.breadth=tot.breadth/hyp.ordi.breadth(dat,grouping=rep(1,dim(dat)[2]),dist.method=dist.method)
		group.vectors<-group.vectors=="YES"
				return(list(species=rownames(dat),eig=ug$eig,tot.breadth=tot.breadth,scaled.breadth=scaled.breadth,distances=distances,group.vectors=group.vectors,centroids.group=centroids.group,plants.ord=ug$vectors,dist.method=dist.method))		
		}
