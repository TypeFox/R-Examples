extract_clinvar <-
function(keyword, localPDB=paste(getwd(), "localPDB",sep="/"), type="both",
        HPO.disease=NULL, genelist=NULL){
    morbidmap=paste(localPDB,"morbidmap",sep="/")
    if(file.exists(localPDB)){
         if(file.exists(paste(localPDB,"variant_summary.txt.gz",sep="/"))){
             clinvar <- paste(localPDB,"variant_summary.txt.gz",sep="/")
             }else{
                 clinvar <- NULL
         }        

         if(file.exists(paste(localPDB,"gene_condition_source_id",sep="/"))){
             gene2dis <- paste(localPDB,"gene_condition_source_id",sep="/")
             }else{
                 gene2dis <- NULL
         }        

         if(file.exists(paste(localPDB,"NAMES.csv.gz",sep="/"))){
             medgene.names <- paste(localPDB,"NAMES.csv.gz",sep="/")
             }else{
                 medgene.names <- NULL
         }        

         if(file.exists(paste(localPDB,"GRtitle_shortname_NBKid.txt",sep="/"))){
             genereview <- paste(localPDB,"GRtitle_shortname_NBKid.txt",sep="/")
             }else{
                 genereview <- NULL
         }        

        }else{
             clinvar <- NULL
             gene2dis <- NULL
             medgene.names <- NULL
             genereview <- NULL
    }     

    #check clinvar database
    if(is.null(clinvar)){
       clinvar <- "ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
       download.path <- paste(getwd(),"localPDB",sep="/")
       if(!file.exists(download.path))
          dir.create(download.path )
       options(timeout = 300)
       if( !file.exists(paste(download.path,"variant_summary.txt.gz",sep="/")))
           download.file(clinvar,
               paste(download.path,"variant_summary.txt.gz",sep="/"),method="auto")
       clinvar <- paste(download.path,"variant_summary.txt.gz",sep="/")
    }

    if(substr(clinvar,nchar(clinvar)-1,nchar(clinvar)) == "gz"){
        clinvar <- read.delim(gzfile(clinvar))
        }else{
            clinvar <- read.delim(clinvar)
    }       

    if(is.null(gene2dis)){
       gene2dis <- "ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/gene_condition_source_id"
       download.path <- paste(getwd(),"localPDB",sep="/")
       if(!file.exists(download.path))
           dir.create(download.path )
       options(timeout = 300)
       if(!file.exists(paste(download.path,"gene_condition_source_id",sep="/")) )
           download.file(gene2dis,paste(download.path,"gene_condition_source_id",sep="/"),method="auto")
       gene2dis <- paste(download.path,"gene_condition_source_id",sep="/")
    }
       gene2dis <- read.delim(gene2dis)
       
    # check medgene 
     if(is.null(medgene.names)){
       medgene.names <- "ftp://ftp.ncbi.nlm.nih.gov/pub/medgen/csv/NAMES.csv.gz"
       download.path <- paste(getwd(),"localPDB",sep="/")
       if(!file.exists(download.path))
          dir.create(download.path )
       options(timeout = 300)
       if( !file.exists(paste(download.path,"NAMES.csv.gz",sep="/")))
           download.file(medgene.names,paste(download.path,"NAMES.csv.gz",sep="/"),method="auto")
       medgene.names <- paste(download.path,"NAMES.csv.gz",sep="/")
    }

    if(substr(medgene.names,nchar(medgene.names)-1,nchar(medgene.names)) == "gz"){
        medgene.names <- read.delim(gzfile(medgene.names))
        }else{
            medgene.names <- read.delim(medgene.names)
    }       

   ## HPO
       if(is.null(HPO.disease)){
          HPO.disease.check <- pheno_extract_HPO(keyword= keyword)
          HPO.disease <- as.character(unique(HPO.disease.check[grep("OMIM",HPO.disease.check[,1]),1]))
       }

  ##check genereview----------------
     if(is.null(genereview)){
       genereview <- "ftp://ftp.ncbi.nlm.nih.gov/pub/GeneReviews/GRtitle_shortname_NBKid.txt"
       download.path <- paste(getwd(),"localPDB",sep="/")
       if(!file.exists(download.path))
          dir.create(download.path )
       options(timeout = 300)
       if( !file.exists(paste(download.path,"GRtitle_shortname_NBKid.txt",sep="/")))
           download.file(genereview,paste(download.path,"GRtitle_shortname_NBKid.txt",sep="/"),method="auto")
       genereview <- paste(download.path,"GRtitle_shortname_NBKid.txt",sep="/")
    }
        genereview <- read.delim(genereview)

   ##function: extract the phenotypeID from the colnames of PhenotypeIDs in clinvar_summary file
   ##databaseID="MedGen","OMIM","GeneReviews",unavailable: "SNOMED CT","Orphanet",
   extract.PhenotypeID <- function(x,database){
      #x=clinvar.extr[49,"PhenotypeIDs"];database="OMIM"
       x.split <- setdiff(unlist(strsplit(as.character(x),";")),"")
       ids <- grep(database,x.split,ignore.case = TRUE)
       if(length(ids)==0){
          id.extr <- ""
          }else if(length(ids)==1){
             x.split.2 <- unlist(strsplit(x.split[ids],","))
             id.extr <- x.split.2[grep(database,x.split.2,ignore.case = TRUE)]
             if(length(id.extr)>1){
                  id.extr <- unique(id.extr)
                  if(length(id.extr)>1){
                      id.extr <- paste(id.extr,collapse=",")
                      }
                  }
             }else{
              id.extr <- c()
              for(j in ids){
                  x.split.2 <- unlist(strsplit(x.split[j],","))
                  id.extr.j <- x.split.2[grep(database,x.split.2,ignore.case = TRUE)]
                  id.extr <- c(id.extr,id.extr.j)
              }  
              id.extr <- paste(id.extr,collapse=",")
       }     
       return(id.extr) 
   }
##########################

 ##begin to search 
       if(!is.null(keyword)){
          gene2dis.d <- gene2dis[grep_split(keyword,gene2dis[,"DiseaseName"]),]
          pheno.yes <- as.character(gene2dis.d[,"DiseaseName"])
          }else if((is.null(keyword)& !is.null(HPO.disease)) | (is.null(keyword) & !is.null(genelist))){
            gene2dis.d <- c()
            pheno.yes <- c()  
            }else{
              gene2dis.d <- gene2dis
              pheno.yes <- c()
       }           
     
       if(!is.null(HPO.disease)){
          HPO.disease.no <- unlist(lapply(HPO.disease,function(x) unlist(strsplit(x,"OMIM:"))[2]))
          gene2dis.d2 <- gene2dis[is.element(gene2dis[,"DiseaseMIM"],HPO.disease.no),]
          pheno.yes2 <- as.character(gene2dis.d2[,"DiseaseName"])
          gene2dis.d <- rbind(gene2dis.d,gene2dis.d2)
          pheno.yes <- union(pheno.yes,pheno.yes2)
       }
   
   #for a given genelist
       if(!is.null(genelist)){
          gene2dis.d3 <- gene2dis[is.element(gene2dis[,"GeneSymbol"],genelist),]
          gene2dis.d <- rbind(gene2dis.d,gene2dis.d3)
       }

       gene2dis.extr <- unique(gene2dis.d)   
      gene2dis.extr$pheno.check <- "no"
      gene2dis.extr[is.element(gene2dis.extr$DiseaseName,pheno.yes),"pheno.check"] <-  "yes"
           
      genes <- unique(as.character(gene2dis.extr[,2]))
     mim.id.pheno <- unique(gene2dis.extr[gene2dis.extr[,"pheno.check"] == "yes",7])
     mim.id.pheno <- mim.id.pheno[!is.na(mim.id.pheno)]
     
      morbidmap <- read.table(morbidmap,header= FALSE,sep="|",quote = "")           
      colnames(morbidmap) <- c("disease","gene","gene.mim.no","location")
     
     ##extract the variants in the genes
     clinvar.extr <- clinvar[is.element(clinvar[,"GeneSymbol"],genes),]
     ##check the phenotype based on the phenotype OMIM ID 
     id.grep <- paste("OMIM",mim.id.pheno,sep=":")
     ##extract mimID
     clinvar.extr$omimID <-  unlist(lapply(clinvar.extr[,"PhenotypeIDs"],function(x) extract.PhenotypeID (x,database="OMIM")))
     clinvar.extr$omim.phenotype <- clinvar.extr$omimID
     phenoid.omim <- unique(clinvar.extr$omimID)
     phenoid.omim <- setdiff(phenoid.omim,c("-",""))
     for(i in phenoid.omim){
         ##i="OMIM:613307"
         if(length(unlist(strsplit(i,","))) == 1){
            id.omim.i <- unlist(strsplit(i,":"))[2]
            disease.i <- as.character(morbidmap[grep(id.omim.i,morbidmap[,1],ignore.case = TRUE),1])
            if(length(disease.i)==0){disease.i=NA}
            }else if(length(unlist(strsplit(i,","))) > 1){
                ids.i <- unlist(strsplit(i,","))
                disease.i <- c()
                for(j in ids.i){
                   ##j=ids.i[1]
                   id.omim.j <- unlist(strsplit(j,":"))[2]
                   disease.j <- as.character(morbidmap[grep(id.omim.j,morbidmap[,1],ignore.case = TRUE),1])
                   disease.i <- c(disease.i,disease.j)
                }
                disease.i <- paste(disease.i,collapse=";")
         }
       if(length( disease.i) >1) {disease.i = paste(disease.i,collapse=";")}
       clinvar.extr[clinvar.extr$omimID==i,"omim.phenotype"] = disease.i
     }

     
     ##extract MedGenID    
     clinvar.extr$medgenID <-  unlist(lapply(clinvar.extr[,"PhenotypeIDs"],function(x) extract.PhenotypeID (x,database="MedGen")))
     clinvar.extr$medgen.phenotype <- clinvar.extr$medgenID
     phenoid.medgen <- unique(clinvar.extr$medgenID)
     phenoid.medgen <- setdiff(phenoid.medgen,"-")
     for(i in phenoid.medgen){
         ##i="medgen:613307"
         if(length(unlist(strsplit(i,","))) == 1){
            id.medgen.i <- unlist(strsplit(i,":"))[2]
            disease.i <- as.character(medgene.names[grep(id.medgen.i,medgene.names[,1],ignore.case = TRUE),2])
            if(length(disease.i)==0){disease.i=NA}
            }else if(length(unlist(strsplit(i,","))) > 1){
                ids.i <- unlist(strsplit(i,","))
                disease.i <- c()
                for(j in ids.i){
                   ##j=ids.i[1]
                   id.medgen.j <- unlist(strsplit(j,":"))[2]
                   disease.j <- as.character(medgene.names[grep(id.medgen.j,medgene.names[,1],ignore.case = TRUE),2])
                   disease.i <- c(disease.i,disease.j)
                }
                disease.i <- paste(disease.i,collapse=";")
         }
       if(length( disease.i) >1) {disease.i <- paste(disease.i,collapse=";")}
       clinvar.extr[clinvar.extr$medgenID==i,"medgen.phenotype"] = disease.i
     }
    
     ##extract GeneReviews
     clinvar.extr$GeneReviewsID <-  unlist(lapply(clinvar.extr[,"PhenotypeIDs"],function(x) extract.PhenotypeID (x,database="GeneReviews")))
     clinvar.extr$GeneReviews.phenotype <- clinvar.extr$GeneReviewsID
     phenoid.GeneReviews <- unique(clinvar.extr$GeneReviewsID)
     phenoid.GeneReviews <- setdiff(phenoid.GeneReviews,"")
     for(i in phenoid.GeneReviews){
         ##i="GeneReviews:NBK1186"
         if(length(unlist(strsplit(i,","))) == 1){
            id.GeneReviews.i <- unlist(strsplit(i,":"))[2]
            disease.i <- unique(as.character(genereview[grep(id.GeneReviews.i,genereview[,3],ignore.case = TRUE),2]))
            if(length(disease.i)==0){disease.i <- NA}
            }else if(length(unlist(strsplit(i,","))) > 1){
                ids.i <- unlist(strsplit(i,","))
                disease.i <- c()
                for(j in ids.i){
                   ##j=ids.i[1]
                   id.GeneReviews.j <- unlist(strsplit(j,":"))[2]
                   disease.j <- unique(as.character(genereview[grep(id.GeneReviews.j,genereview[,3],ignore.case = TRUE),2]))
                   disease.i <- c(disease.i,disease.j)
                }
                disease.i <- paste(disease.i,collapse=";")
         }
       if(length( disease.i) >1) {disease.i <- paste(disease.i,collapse=";")}
       clinvar.extr[clinvar.extr$GeneReviewsID==i,"GeneReviews.phenotype"] = disease.i
     }
     extract <- list(gene2dis.extr,clinvar.extr)
     names(extract) <- c("gene2dis","variants")
     return(extract)
}
