

identifyLncGraphW<-function(moleculeList,graphList,type="gene_lncRNA",background=getBackgroundLnc(type),
   order="pvalue",decreasing=FALSE,bet=1){
     #require(BiasedUrn)
	 lnc2Name<-get("lnc2Name",envir=envData)
	 if(typeof(moleculeList)!="character"&&!is.null(moleculeList)){
		print("warning: your moleculeList must be 'character' vector. Because the type of your current moleculeList is not correct, it has been conveted arbitrarily using the function as.character().")
        moleculeList<-as.character(unlist(moleculeList))
	  }
   ############################################################
      if(!exists("envData")) envData<-initialize()
	  graphList.length<-length(graphList)
	  if(graphList.length==0){
	     print("warning: there is no any graphs in graphList or these graphs are not available for pathway analyses.")
	  }	  
	 
	     if(graphList.length>0){
	         gene2path<-get("gene2path",envir=envData)
	         org.path<-unique(as.character(gene2path[,2]))
		     org.path<-substring(org.path,nchar(org.path)-4,nchar(org.path))
		     graphList<-graphList[sapply(graphList,function(x) substring(x$number,0,5) %in% org.path)]
		 }
	
      annList<-list()
	  if(graphList.length>0){
      for(i in 1:length(graphList)){
            ann<-list(pathwayId=character(),pathwayName="not known",annMoleculeList=character(),annMoleculeNumber=0,
                      annBgMoleculeList=character(),annBgNumber=0,moleculeNumber=0,weight=1,bgNumber=0,pvalue=1,fdr=1)

			ann$pathwayId<-paste("path:",graphList[[i]]$number,sep="")
            KOList<-unique(unlist(strsplit(unlist(strsplit(V(graphList[[i]])$names," ")),";")))
            
			#KOList<-unique(unlist(strsplit(V(graphList[[i]])$names,"[ ;]")))
   
        ###########################start###################################
        if(type=="gene"||type=="gene_lncRNA"){	
		#if(length(grep("[A-Z]",KOList)))
            graphGeneList <- unique(KOList[-grep("[A-Z]",KOList)])
			#else graphGeneList<-KOList
       }	
	   
	   if(type=="lncRNA"||type=="gene_lncRNA"){	
            graphMirnaList<-unique(KOList[grep("[A-Z]",KOList)]) 
            graphMirnaList<-intersect(as.character(graphMirnaList),as.character(lnc2Name[,2]))
	   }
	  
	   if(type=="gene_lncRNA"){
	        graphMoleculeList<-unique(c(graphGeneList,graphMirnaList))
	   }
	   else if(type=="gene"){
	        graphMoleculeList<-graphGeneList   
	   }
	   
	   else if(type=="lncRNA"){
	        graphMoleculeList<-graphMirnaList  	   
	   }
	   ##############################end ##################################
	   background<-getBackgroundLnc("gene")
       annotatedMoleculeList<-intersect(graphMoleculeList,moleculeList)
	   annotatedBackgroundList<-intersect(graphMoleculeList,getBackgroundLnc(type))	
	   if(length(annotatedMoleculeList[grep("[A-Z]",annotatedMoleculeList)])){ 
	   GeneNum<-length(graphGeneList)
	   edgelist<-get.edgelist(graphList[[i]])
	   if(type=="lncRNA"||type=="gene_lncRNA"){
	       hit<-c(edgelist[which(edgelist[,1]%in%intersect(graphMirnaList,moleculeList)),2],edgelist[which(edgelist[,2]%in%intersect(graphMirnaList,moleculeList)),1])
	       LncLinkGeneNum<-length(unique(hit))
	   }
            
            pathwayName<-graphList[[i]]$title
            if(length(pathwayName)!=0)
                ann$pathwayName<-pathwayName
            ann$annMoleculeList<-intersect(graphGeneList,moleculeList) 
         
            ann$annMoleculeNumber<-length(ann$annMoleculeList)
			ann$moleculeNumber<-length(intersect(moleculeList,getBackgroundLnc("gene")))
			ann$annBgMoleculeList<-intersect(graphMoleculeList,getBackgroundLnc(type)) 
            ann$annBgNumber<-length(intersect(graphGeneList,getBackgroundLnc("gene")))

            ann$bgNumber<-length(getBackgroundLnc("gene"))
			if(LncLinkGeneNum){
            ann$weight<-1+bet*(-log2(LncLinkGeneNum/GeneNum))
			#ann$weight<-1
			error1<-try( ann$pvalue<-pWNCHypergeo(ann$annMoleculeNumber-1,ann$annBgNumber,
                 ann$bgNumber-ann$annBgNumber,ann$moleculeNumber,ann$weight,precision=1E-100,lower.tail=FALSE),TRUE)
              if(!(error1<=1)){
                 print(graphList[[i]]$number)
				 print( ann$pvalue<-1-phyper(ann$annMoleculeNumber-1,ann$annBgNumber,
                 ann$bgNumber-ann$annBgNumber,ann$moleculeNumber))
				 ann$pvalue<-1-phyper(ann$annMoleculeNumber-1,ann$annBgNumber,
                 ann$bgNumber-ann$annBgNumber,ann$moleculeNumber)
              }	

			}
		
			
			#print(LncLinkGeneNum)
			#print(GeneNum)
			#print(ann$weight)
			
            annList[[i]]<-ann
			}
			
			else{ 
			 pathwayName<-graphList[[i]]$title
            if(length(pathwayName)!=0)
                ann$pathwayName<-pathwayName
            ann$annMoleculeList<-intersect(graphGeneList,moleculeList) 
         
            ann$annMoleculeNumber<-length(ann$annMoleculeList)
			ann$moleculeNumber<-length(intersect(moleculeList,getBackgroundLnc("gene")))
			ann$annBgMoleculeList<-intersect(graphMoleculeList,getBackgroundLnc(type)) 
            ann$annBgNumber<-length(intersect(graphGeneList,getBackgroundLnc("gene")))

            ann$bgNumber<-length(getBackgroundLnc("gene"))
			
			}
	     annList[[i]]<-ann
      } 
	  }
	  fdr.est<-function(p)
{
    m <- length(ind <- which(!is.na(p)))
    fdr <- rep(NA, length(p))
    stat <- cbind(1:length(p), p, fdr)
    stat[ind, 3] <- unlist(lapply(stat[ind, 2], function(x) {
        c <- length(which(stat[ind, 2] <= x))
        m * x/c
    }))
    stat[ind, ] <- stat[ind, ][order(stat[, 2], decreasing = TRUE), 
        ]
    stat[ind, 3] <- cummin(stat[ind, 3])
    fdr <- stat[order(stat[, 1]), 3]
    fdr[which(fdr > 1)] <- 1
    return(fdr)
}
	  if(length(annList)>0){
	     p_value<-sapply(annList,function(x) return(x$pvalue))
         #fdrtool.List<-fdrtool(p_value,statistic="pvalue",plot=FALSE,verbose=FALSE)
         	 
         #print(fdrtool.List$qval)
         #for(i in seq(annList)){
         #   annList[[i]]$qvalue<-fdrtool.List$qval[i]
		 #	annList[[i]]$lfdr<-fdrtool.List$lfdr[i]
         #}
		 fdr.List<-fdr.est(p_value)
		 for(i in seq(annList)){
		     annList[[i]]$fdr<-fdr.List[i]
		 }
         #names(annList)<-sapply(graphList,function(x) x$number)
         annList<-annList[sapply(annList,function(x) x$annMoleculeNumber>0)]
         annList<-annList[order(sapply(annList,function(x) x[[order]]),decreasing=decreasing)]   
	  }
	  return(annList)	

}

