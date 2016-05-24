extract_genes_orphanet <-
function(keyword, localPDB = paste(getwd(),"localPDB",sep="/"), HPO.disease = NULL){
    if(file.exists(localPDB)){
         if(file.exists(paste(localPDB,"en_product6.xml",sep="/"))){
             orphanet <- paste(localPDB,"en_product6.xml",sep="/")
             }else{
                 orphanet <- NULL
         }        
        }else{
             orphanet <- NULL 
    }     
   
   # HPO
       if(is.null(HPO.disease)){
          HPO.disease.check <- pheno_extract_HPO(keyword= keyword)
          HPO.disease <- as.character(unique(HPO.disease.check[grep("ORPHANET",HPO.disease.check[,1]),1]))
       }
 
    if(is.null(orphanet)){
       orphanet <- "http://www.orphadata.org/data/xml/en_product6.xml"
       download.path = paste(getwd(),"localPDB",sep="/")
       if(!file.exists(download.path))
           dir.create(download.path )
       options(timeout = 300)
       if(!file.exists(paste(download.path,"en_product6.xml",sep="/")))
            download.file(orphanet,paste(download.path,"en_product6.xml",sep="/"),method="auto")
       orphanet <- paste(download.path,"en_product6.xml",sep="/")
    }

    doc <- xmlTreeParse(orphanet,useInternalNodes = TRUE)        
    nodes <- getNodeSet(doc, "//Disorder")
    lists <- nodesToList(nodes)
    
    pheno2gene <- c()
    for(i in 1:length(lists)){     
       OrphaNumber.i <- lists[i][[1]]$OrphaNumber
       Phenotype.i <- lists[i][[1]]$Name$text
       if(is.element(paste("ORPHANET",OrphaNumber.i,sep=":"),HPO.disease) | length(grep_split(keyword,Phenotype.i))>0){
          for(j in 1:(length(lists[i][[1]]$GeneList)-1)){
              GeneType.i.j <- lists[i][[1]]$GeneList[j]$Gene$GeneType$Name$text
              GeneName.i.j <- lists[i][[1]]$DisorderGeneAssociationList[j]$DisorderGeneAssociation$Gene$Name$text
              GeneSymbol.i.j <- lists[i][[1]]$DisorderGeneAssociationList[j]$DisorderGeneAssociation$Gene$Symbol
              AssociationType.i.j <- lists[i][[1]]$DisorderGeneAssociationList[j]$DisorderGeneAssociation$DisorderGeneAssociationType$Name$text
              AssociationStatus.i.j <- lists[i][[1]]$DisorderGeneAssociationList[j]$DisorderGeneAssociation$DisorderGeneAssociationStatus$Name$text
              pheno2gene <- rbind(pheno2gene,c(OrphaNumber.i,Phenotype.i,GeneSymbol.i.j,GeneName.i.j,GeneType.i.j,AssociationType.i.j,AssociationStatus.i.j))
          }
       }   
    }
    colnames(pheno2gene) <- c("OrphaNumber","Phenotype","GeneSymbol","GeneName","GeneType","AssociationType","AssociationStatus")  
    return(pheno2gene)
}
