
getBackgroundLnc<-function(type="gene_lncRNA"){
      if(!exists("envData")) envData<-initialize()
	  if(type=="gene"||type=="gene_lncRNA"){
	  keggGene2gene<-get("keggGene2gene",envir=envData) 
	  background<-unique(as.character(keggGene2gene[,2]))
	  geneBackground<-sapply(strsplit(background,":"),function(x) x[2])
	  newBackground<-geneBackground
	  }
	  if(type=="lncRNA"||type=="gene_lncRNA"){
	  lncBackground<-get("lncBackground",envir=envData)
      newBackground<-lncBackground
	  }
	  if(type=="gene_lncRNA"){
	  newBackground<-union(geneBackground,lncBackground)
	  }
	 return(newBackground)
}