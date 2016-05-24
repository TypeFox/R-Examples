pheno_extract_HPO <-
function(keyword, localPDB = paste(getwd(),"localPDB",sep="/")){
    if(file.exists(localPDB)){
         if(file.exists(paste(localPDB,"phenotype_annotation.tab",sep="/"))){
             HPO <- paste(localPDB,"phenotype_annotation.tab",sep="/")
             }else{
                 HPO <- NULL
         }        
         if( file.exists(paste(localPDB,"diseases_to_genes.txt",sep="/"))) {  
             diseases_to_genes <- paste(localPDB,"diseases_to_genes.txt",sep="/")
             }else{
                  diseases_to_genes <- NULL
         }         
        }else{
            HPO <- NULL; diseases_to_genes <- NULL  
    }     
         
    #check hgnc database
    if(is.null(HPO)){
       HPO <- "http://compbio.charite.de/hudson/job/hpo.annotations/lastStableBuild/artifact/misc/phenotype_annotation.tab"
       download.path <- paste(getwd(),"localPDB",sep="/")
       if(!file.exists(download.path))
          dir.create(download.path )
       options(timeout = 300)
       if( !file.exists(paste(download.path,"phenotype_annotation.tab",sep="/")))
           download.file(HPO,paste(download.path,"phenotype_annotation.tab",sep="/"),method="auto")
       HPO <- paste(download.path,"phenotype_annotation.tab",sep="/")
    }
     HPO <- read.delim(HPO,header= FALSE)
     
     #input disease2gene dataset
    if(is.null(diseases_to_genes)){
       diseases_to_genes <- "http://compbio.charite.de/hudson/job/hpo.annotations.monthly/lastStableBuild/artifact/annotation/diseases_to_genes.txt"
       download.path <- paste(getwd(),"localPDB",sep="/")
       if(!file.exists(download.path))
          dir.create(download.path )
       options(timeout = 300)
       if( !file.exists(paste(download.path,"diseases_to_genes.txt",sep="/")))
           download.file(diseases_to_genes,paste(download.path,"diseases_to_genes.txt",sep="/"),method="curl")
       diseases_to_genes <- paste(download.path,"diseases_to_genes.txt",sep="/")
    }
     diseases_to_genes <- read.delim(diseases_to_genes,fill= TRUE , flush= TRUE,header=FALSE,col.names=c("diseaseId","geneID","GeneSymbol"),comment.char = "#")
    
    if(!is.null(keyword)){ 
      HPO.merge <- unlist(apply(HPO,1,function(x) paste(as.character(x),collapse="_")))
      HPO.j <- HPO[grep_split(keyword,HPO.merge),]

      HPO.j.sim <-  unique(HPO.j[,c(1:3,6,12)])
      colnames(HPO.j.sim) <- c("Database","ID","DiseaseName","reference","Synonym")
          
      db.id <- paste(HPO.j.sim[,1],HPO.j.sim[,2],sep=":")
      diseases_to_genes.j <- diseases_to_genes[is.element(diseases_to_genes[,1],db.id),]
      colnames(diseases_to_genes.j) <- c("DiseaseID","GeneID","GeneName")
      diseases_to_genes.j$Synonym <- diseases_to_genes.j$DiseaseName <-  ""
      for(i in unique(diseases_to_genes.j[,1])){
          diseases_to_genes.j[diseases_to_genes.j[,1] == i,"DiseaseName"] <- as.character(unique(HPO.j.sim[db.id==i,3]))
          diseases_to_genes.j[diseases_to_genes.j[,1] == i,"Synonym"] <- as.character(unique(HPO.j.sim[db.id==i,5]))         
      }
       
      return(diseases_to_genes.j)  
      }else{
      return(list(HPO,diseases_to_genes))
    } 
}
