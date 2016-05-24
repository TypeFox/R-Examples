genes_compile <-
function(HPO,orphanet,omim,clinvar,uniprot,localPDB = paste(getwd(),"localPDB",sep="/")){

# trim the gene names: 'ORF' -> 'orf'
    gene.orf <- function(x){
       # x = "C8ORF37"
        if(length(grep("ORF",x))>0){
          x <- paste(unlist(strsplit(x,"ORF")),collapse="orf")
          }
        return(x)    
    }
    
    hgnc <- read.delim(gzfile(paste(localPDB,"hgnc_complete_set.txt.gz",sep="/")))
    refFlat <- read.delim(gzfile(paste(localPDB,"refFlat.txt.gz",sep="/")),header= FALSE)    
    genes <- unique(c(as.character(HPO[,3]), as.character(orphanet[,3]), as.character(omim[,6]), as.character(clinvar[,2]), as.character(uniprot[,1])))
    genes <- genes[genes != "" & genes != "missing"]
    genes.trim <- unique(unlist(lapply(genes,gene.orf)))
    refFlat.extract <- refFlat[is.element(refFlat[,1],genes.trim),]
    gene.position <- matrix(,nrow=length(genes.trim),ncol=5)
    colnames(gene.position) <- c("Gene.symbol","chr","strand","start","end")
    rownames(gene.position) <- gene.position[,1] <- genes.trim
    for(g in genes.trim){
       refFlat.extract.g <- refFlat.extract[refFlat.extract[,1] == g,]
       if(nrow(refFlat.extract.g) > 0 ){
           gene.position[g,"chr"] <- paste(unique(refFlat.extract.g[,3]),collapse=",")
           gene.position[g,"strand"] <- paste(unique(refFlat.extract.g[,4]),collapse=",")
           gene.position[g,"start"] <- min(refFlat.extract.g[,5])
           gene.position[g,"end"] <- max(refFlat.extract.g[,6])      
       }
    }
    
    HPO[,3] <- unlist(lapply(as.character(HPO[,3]),gene.orf))
    orphanet[,3] <- unlist(lapply(as.character(orphanet[,3]),gene.orf))
    omim[,6] <- unlist(lapply(as.character(omim[,6]),gene.orf))
    clinvar[,2] <- unlist(lapply(as.character(clinvar[,2]),gene.orf))
    uniprot[,1] <- unlist(lapply(as.character(uniprot[,1]),gene.orf))
    gene2pheno <- matrix(,nrow=length(genes.trim),ncol=8)
    colnames(gene2pheno) <- c("Entrez.Gene.ID","Approved.Name","Synonyms","HPO","Orphanet","OMIM","ClinVar","Uniprot")
    rownames(gene2pheno) <- genes.trim
    hgnc.extract <- hgnc[is.element(hgnc$Approved.Symbol,genes.trim),]
    rownames(hgnc.extract) <- hgnc.extract$Approved.Symbol
    gene2pheno[,c("Synonyms")] <- as.character(hgnc.extract[genes.trim,c("Synonyms")])
    gene2pheno[,c("Approved.Name")] <- as.character(hgnc.extract[genes.trim,c("Approved.Name")])
    gene2pheno[,"Entrez.Gene.ID"] <- as.character(hgnc.extract[genes.trim,c("Entrez.Gene.ID")])       
      for(i in genes.trim){
       # i = genes[1]
       gene2pheno[i,"HPO"] <- paste(unique(HPO[HPO[,3] == i,4]),collapse=";")
       gene2pheno[i,"Orphanet"] <- paste(orphanet[orphanet[,3] == i,2],collapse=";")
       if(!is.null(omim))
           gene2pheno[i,"OMIM"] <- paste(omim[omim[,6] == i,1],collapse=";")
           gene2pheno[i,"ClinVar"] <- paste(clinvar[clinvar[,2] == i,4],collapse=";")
           gene2pheno[i,"Uniprot"] <- paste(uniprot[uniprot[,1] == i,2],collapse=";")
    }
    gene2pheno <- cbind(gene.position,gene2pheno)
    
    #add the omim missing genes
    if(!is.null(omim)){
       omim.missing <- omim[omim$Approved.Symbol == "missing",]
       if(nrow(omim.missing) > 0){
         gene2pheno <- rbind(gene2pheno,matrix(,nrow=nrow(omim.missing),ncol=8+5))
         gene2pheno[(nrow(gene2pheno)-nrow(omim.missing)+1):nrow(gene2pheno),"Synonyms"] <- as.character(omim.missing[,"gene"])
         gene2pheno[(nrow(gene2pheno)-nrow(omim.missing)+1):nrow(gene2pheno),"OMIM"] <- as.character(omim.missing[,"disease"])          
      }
    }
    return(gene2pheno)
}

