getSymbolFromGene<-function(geneList){
getKGeneFromGene<-function(geneList){
	  geneList<-as.character(geneList)
      if(!exists("envData")) envData<-initialize()
	  keggGene2gene<-get("keggGene2gene",envir=envData)
keggGeneList<-unique(as.character(keggGene2gene[as.character(keggGene2gene[,2]) %in% paste("ncbi-geneid",geneList,sep=":"),1]))
      return(keggGeneList)
}
getSymbolFromKGene<-function(keggGeneList){
	  keggGeneList<-as.character(keggGeneList)
      if(!exists("envData")) envData<-initialize()
      gene2symbol<-get("gene2symbol",envir=envData)
      symbolList<-unique(as.character(sapply(strsplit(as.character(gene2symbol[as.character(gene2symbol[,1]) %in% keggGeneList,2]),":"),function(x) return (x[2]))))
      return(symbolList)
}
      return(getSymbolFromKGene(getKGeneFromGene(geneList)))
}