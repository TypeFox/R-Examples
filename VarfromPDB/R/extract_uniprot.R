extract_uniprot <-
function(keyword, localPDB = paste(getwd(),"localPDB",sep="/"), HPO.disease = NULL, genelist = NULL){
    if(file.exists(localPDB)){
         if(file.exists(paste(localPDB,"humsavar.txt",sep="/"))){
             uniprot <- paste(localPDB,"humsavar.txt",sep="/")
             }else{
                 uniprot <- NULL
         }        
        }else{
             uniprot <- NULL
    }     

   # HPO
       if(is.null(HPO.disease)){
          HPO.disease.check <- pheno_extract_HPO(keyword= keyword)
          HPO.disease <- as.character(unique(HPO.disease.check[grep("OMIM",HPO.disease.check[,1]),1]))
       }
 
    if(is.null(uniprot)){
       uniprot <- "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/variants/humsavar.txt"
       download.path <- paste(getwd(),"localPDB",sep="/")
       if(!file.exists(download.path))
          dir.create(download.path )
       options(timeout = 300)
       if( !file.exists(paste(download.path,"humsavar.txt",sep="/")))
           download.file(uniprot,paste(download.path,"humsavar.txt",sep="/"),method="auto")
       uniprot <- paste(download.path,"humsavar.txt",sep="/")
    }

    dat <- readLines(uniprot,n=300)
    num.1 <- lapply(dat,function(x) grep("Swiss-Prot",x))
    num.2 <- lapply(dat,function(x) grep("dbSNP",x))
    names(num.1) <- names(num.2) <- 1:length(dat)
    skip.num <- intersect(as.numeric(names(unlist(num.1)))+2,as.numeric(names(unlist(num.2)))+1)
    dat <- read.fwf(uniprot, widths = c(10,11,12,15,14,12,200),skip=skip.num,n = -1,row.names=NULL)  
    colnames(dat) <- c("GeneSymbol","swiss.prot.AC","FTId","AA.change","type","dbSNP","DiseaseName")
    dat$DiseaseMIM <- unlist(lapply(dat$DiseaseName,function(x) unlist(strsplit(unlist(strsplit(as.character(x),"MIM:"))[2],"]"))))
    dat[,"GeneSymbol"] <- unlist(lapply(dat[,"GeneSymbol"],function(x) str_trim(as.character(x), side = c("both"))))
    
#extract
       if(!is.null(keyword)){
          dat.d <- dat[grep_split(keyword,dat[,"DiseaseName"]),]
          pheno.yes <- as.character(dat.d[,"DiseaseName"])
          }else if((is.null(keyword)& !is.null(HPO.disease)) | (is.null(keyword) & !is.null(genelist))){
            dat.d <- c()
            pheno.yes <- c()  
            }else{
              dat.d <- dat
              pheno.yes <- c()
       }           

       if(!is.null(HPO.disease)){
          HPO.disease.no <- unlist(lapply(HPO.disease,function(x) unlist(strsplit(x,"OMIM:"))[2]))
          dat.d2 <- dat[is.element(dat[,"DiseaseMIM"],HPO.disease.no),]
          pheno.yes2 <- as.character(dat.d2[,"DiseaseName"])
          dat.d <- rbind(dat.d,dat.d2)
          pheno.yes <- union(pheno.yes,pheno.yes2)
       }

       if(!is.null(genelist)){
          dat.d3 <- dat[is.element(as.character(dat[,"GeneSymbol"]),genelist),]
          dat.d <- rbind(dat.d,dat.d3)
       }

      dat.extr <- unique(dat.d)   
      dat.extr$pheno.check <- "no"
      dat.extr[is.element(dat.extr$DiseaseName,pheno.yes),"pheno.check"] <-  "yes"      
      genes.extr <- unique(dat.extr[,c("GeneSymbol","DiseaseName","pheno.check")])
      genes.extr <- genes.extr[genes.extr$DiseaseName != "-",]      
      return(list(genes.extr,dat.extr))
}
