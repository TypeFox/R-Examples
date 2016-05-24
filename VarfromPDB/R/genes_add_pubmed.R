genes_add_pubmed <-
function(genepdb,pubmed,localPDB = paste(getwd(),"localPDB",sep="/")){
## filter the pubmed
   pubmed <- pubmed[!is.na(pubmed[,1])&(pubmed[,2]!=""),]

   pubmed_genes <- unlist(lapply(as.character(pubmed[,"Approved.Symbol"]),function(x) unlist(strsplit(x,", "))))
   pdb_genes <- as.character(genepdb[!is.na(genepdb[,"chr"]),"Gene.symbol"])

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
    genes <- unique(c(as.character(pdb_genes), as.character(pubmed_genes)))
    genes <- genes[genes != ""]
    genes <- genes[-grep("missing",genes)]
    genes.trim <- unique(unlist(lapply(genes,gene.orf)))
    genes_add <- setdiff(genes.trim,pdb_genes)
    refFlat.extract <- refFlat[is.element(refFlat[,1],genes.trim),]
    gene.position <- matrix(,nrow=length(genes_add),ncol=8)
    colnames(gene.position) <- c("Gene.symbol","chr","strand","start","end","Entrez.Gene.ID","Approved.Name","Synonyms")
    rownames(gene.position) <- gene.position[,1] <- genes_add
    for(g in genes_add){
       refFlat.extract.g <- refFlat.extract[refFlat.extract[,1] == g,]
       if(nrow(refFlat.extract.g) > 0 ){
           gene.position[g,"chr"] <- paste(unique(refFlat.extract.g[,3]),collapse=",")
           gene.position[g,"strand"] <- paste(unique(refFlat.extract.g[,4]),collapse=",")
           gene.position[g,"start"] <- min(refFlat.extract.g[,5])
           gene.position[g,"end"] <- max(refFlat.extract.g[,6])                 
       }
    }
    hgnc.extract <- hgnc[is.element(hgnc$Approved.Symbol,genes_add),]
    rownames(hgnc.extract) <- hgnc.extract$Approved.Symbol
    gene.position[rownames(hgnc.extract),c("Synonyms")] <- as.character(hgnc.extract[rownames(hgnc.extract),c("Synonyms")])
    gene.position[rownames(hgnc.extract),c("Approved.Name")] <- as.character(hgnc.extract[rownames(hgnc.extract),c("Approved.Name")])
    gene.position[rownames(hgnc.extract),c("Entrez.Gene.ID")] <- as.character(hgnc.extract[rownames(hgnc.extract),c("Entrez.Gene.ID")])
   
    geneAll <- matrix(,nrow=nrow(genepdb)+nrow(gene.position),ncol=ncol(genepdb))
    colnames(geneAll) <- colnames(genepdb)
    geneAll[1:nrow(genepdb),] <- as.matrix(genepdb[,])
    geneAll[(1+nrow(genepdb)):nrow(geneAll),1:ncol(gene.position)] <- as.matrix(gene.position) 
    
    check.gene <- function(gene,x){
        x.split <- unlist(strsplit(as.character(x),", "))
        check.i <- intersect(gene,x.split)
        check.result <- ifelse(length(check.i)>0,T,F)
        return(check.result)
    }
    
    geneAll <- cbind(geneAll,"")
    colnames(geneAll)[ncol(geneAll)] <- "pubmed"
    for(i in 1:nrow(geneAll)){
       gene.i <- geneAll[i,1]  
       pheno.i <- as.character(unique(pubmed[unlist(lapply(pubmed[,"Approved.Symbol"],function(x) check.gene(gene.i,x))), "Phenotype"]))
       if(length(pheno.i) > 0){
           geneAll[i,"pubmed"] <- paste(pheno.i,collapse=";")
           }else{
              geneAll[i,"pubmed"] <- ""
       }       
    }    
    return(geneAll)
}

