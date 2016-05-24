getKGeneFromGene<-function(geneList){
	  geneList<-as.character(geneList)
      if(!exists("envData")) envData<-initialize()
	  keggGene2gene<-get("keggGene2gene",envir=envData)
keggGeneList<-unique(as.character(keggGene2gene[as.character(keggGene2gene[,2]) %in% paste("ncbi-geneid",geneList,sep=":"),1]))
      return(keggGeneList)
}