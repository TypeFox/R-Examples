getBackground <-
function(type="gene_miRNA"){
      if(!exists("k2ri")) k2ri<-initializeK2ri()
	  if(type=="gene"||type=="gene_miRNA"){
	  geneBackground<-GetK2riData("BGGene") 
	  newBackground<-geneBackground
	  }
	
	  if(type=="miRNA"||type=="gene_miRNA"){
	  miRNABackground<-GetK2riData("BGMiRNA")
	  newBackground<-miRNABackground
	  }

	  if(type=="gene_miRNA"){
	  newBackground<-union(geneBackground,miRNABackground)
	  }
	  
	 return(newBackground)
}
