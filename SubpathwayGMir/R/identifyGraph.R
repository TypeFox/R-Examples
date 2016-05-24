identifyGraph <-
function(moleculeList,graphList,type="gene_miRNA",background=getBackground(type),order="pvalue",decreasing=FALSE)
{
	 if(typeof(moleculeList)!="character"&&!is.null(moleculeList)){
		print("warning: your moleculeList must be 'character' vector.")
        moleculeList <- as.character(unlist(moleculeList))
	  }
   ############################################################
      if(!exists("k2ri")) k2ri <- initializeK2ri()
	  graphList.length <- length(graphList)
	  if(graphList.length==0){
	     print("warning: there is no any graphs in graphList.")
	  }	  
	  
	  if(graphList.length>0){
	         gene2path<-GetK2riData("gene2path")
	         org.path<-unique(as.character(gene2path[,2]))
		     org.path<-substring(org.path,nchar(org.path)-4,nchar(org.path))
		     graphList<-graphList[sapply(graphList,function(x) substring(x$number,0,5) %in% org.path)]
		 }
	 
      annList<-list()
	  if(graphList.length>0){
      for(i in 1:length(graphList)){
            ann<-list(pathwayId=character(),pathwayName="not known",annMoleculeList=character(),annMoleculeNumber=0,
                      annBgMoleculeList=character(),annBgNumber=0,moleculeNumber=0,bgNumber=0,pvalue=1,fdr=1)

			ann$pathwayId<-paste("path:",graphList[[i]]$number,sep="")
            KOList<-unique(unlist(strsplit(unlist(strsplit(V(graphList[[i]])$names," ")),";")))
            
			#KOList<-unique(unlist(strsplit(V(graphList[[i]])$names,"[ ;]")))
   
        ###########################start###################################
       if(type=="gene"||type=="gene_miRNA"){	
            graphGeneList <- KOList[-grep("-",KOList)]
       }	
       
	   if(type=="miRNA"||type=="gene_miRNA"){	
            graphMirnaList<-KOList[grep("-",KOList)]#graphMirnaList<-KOList[grep("[A-Z]",KOList)] 
	   }
	
	   if(type=="gene_miRNA"){
	        graphMoleculeList<-c(graphGeneList,graphMirnaList)
	   }
	   else if(type=="gene"){
	        graphMoleculeList<-graphGeneList   
	   }
	
	   else if(type=="miRNA"){
	        graphMoleculeList<-graphMirnaList  	   
	   }
	   ##############################end ##################################
	   background<-getBackground(type)
       annotatedMoleculeList<-intersect(graphMoleculeList,moleculeList)
	   annotatedBackgroundList<-intersect(graphMoleculeList,background)	
            
            pathwayName<-graphList[[i]]$title
            if(length(pathwayName)!=0)
                ann$pathwayName<-pathwayName
            ann$annMoleculeList<-annotatedMoleculeList 
         
            ann$annMoleculeNumber<-length(annotatedMoleculeList)
			ann$annBgMoleculeList<-annotatedBackgroundList
            ann$annBgNumber<-length(annotatedBackgroundList)

            ann$moleculeNumber<-length(moleculeList)
            ann$bgNumber<-length(background)

            ann$pvalue<-1-phyper(ann$annMoleculeNumber-1,ann$annBgNumber,
                 ann$bgNumber-ann$annBgNumber,ann$moleculeNumber)
            
            annList[[i]]<-ann
      } 
	  }
	  if(length(annList)>0){
	     p_value<-sapply(annList,function(x) return(x$pvalue))
		 fdr.List<-fdr.est(p_value)
		 for(i in seq(annList)){
		     annList[[i]]$fdr<-fdr.List[i]
		 }
         annList<-annList[sapply(annList,function(x) x$annMoleculeNumber>0)]
         annList<-annList[order(sapply(annList,function(x) x[[order]]),decreasing=decreasing)]   
	  }
	  return(annList)	

}
