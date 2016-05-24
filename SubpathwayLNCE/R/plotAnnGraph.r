plotAnnGraph<-function(pathwayId,graphList,ann,gotoKEGG=FALSE,orgSpecific=TRUE,displayInR=TRUE,match=TRUE,vertex.frame.color="red",...){
     url.list<-c()
	 warning_result<-FALSE
     newPathwayId<-sapply(pathwayId,function(x) unlist(strsplit(x,":"))[2])
     entireNewPathwayId<-substring(newPathwayId,0,5)

	 for(i in seq(pathwayId)){
         if(match==FALSE){
             matchGraph<-graphList[[entireNewPathwayId[i]]]
        }else{
             matchGraph<-graphList[[newPathwayId[i]]]
        }
         org_idType<-unlist(strsplit(matchGraph$org,";"))
         org<-org_idType[1]
		 
         #ann[sapply(ann,function(x) ifelse(x$pathwayId==pathwayId,TRUE,FALSE))]
         annMoleculeList<-ann[sapply(ann,function(x) ifelse(x$pathwayId==pathwayId[i],TRUE,FALSE))][[1]]$annMoleculeList

         
		 KOList<-""
         if(org=="hsa"){
				 
					 KOList<-annMoleculeList
		}
		else{stop(paste("graph ",i,"  error: it is not ec, ko, or current org graph.",sep=""))}	 
		 
		 
		 
		 
		if(displayInR==TRUE){
             componentList<-KOList
             frame.color<-sapply(V(matchGraph)$names, function(x)
                 ifelse(length(intersect(unlist(strsplit(unlist(strsplit(x," ")),";")),
                 componentList))>0,vertex.frame.color,"dimgray"))
             plotGraphL(matchGraph,vertex.frame.color=frame.color,...)
		 }
         if(gotoKEGG==TRUE){
		     limit.length<-250
             if(orgSpecific==TRUE){
	             org<-"hsa"
	             annGeneList<-sapply(strsplit(getKGeneFromGene(annMoleculeList),":"), function(x) x[2])
	             temp<- paste(c(paste(org,entireNewPathwayId[i],sep=""),annGeneList),sep="",collapse="+")
                 url <- paste("http://www.genome.ad.jp/dbget-bin/show_pathway?",temp,sep="")
	             #print(url)
				 if(length(unlist(strsplit(url,"")))>limit.length){
				    warning_result<-TRUE
					url.list<-c(url.list,url)
				    url<-substring(url,0,limit.length)
					print(paste("warning: gene numbers in ", pathwayId[i], " are too large.",sep=""))
				 }
                 browseURL(url)
            }else{	  
	             temp<- paste(c(paste(org,entireNewPathwayId[i],sep=""),KOList),sep="",collapse="+")
                 url <- paste("http://www.genome.ad.jp/dbget-bin/show_pathway?",temp,sep = "")
				 if(length(unlist(strsplit(url,"")))>limit.length){
				    url.list<-c(url.list,url)
					url<-substring(url,0,limit.length)
					print(paste("warning: gene numbers in ", pathwayId[i], " are too large.",sep=""))
				 }
                 browseURL(url)
	        }
        }
	}
	
	if(warning_result==TRUE){
	return(url.list)
	}
	V(matchGraph)$type[is.na((V(matchGraph)$type))]<-"gene"
	return(matchGraph)
}
