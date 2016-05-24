extract_omim <-
function(keyword, omim.apiKey,
         localPDB = paste(getwd(),"localPDB",sep="/"), type = "both", HPO.disease = NULL, genelist = NULL){
    if(file.exists(localPDB)){
         if(file.exists(paste(localPDB,"hgnc_complete_set.txt.gz",sep="/"))){
             hgnc <- paste(localPDB,"hgnc_complete_set.txt.gz",sep="/")
             }else{
                 hgnc <- NULL
         }        
        }else{
             hgnc <- NULL
    }     
       
       morbidmap = paste(localPDB,"morbidmap",sep="/")
       morbidmap <- read.table(morbidmap,header= FALSE,sep="|",quote = "")           
       colnames(morbidmap) <- c("disease","gene","gene.mim.no","location")
    #check hgnc database
    if(is.null(hgnc)){
       hgnc <- "ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc_complete_set.txt.gz"
       download.path <- paste(getwd(),"localPDB",sep="/")
       if(!file.exists(download.path))
          dir.create(download.path )
       options(timeout = 300)
       if( !file.exists(paste(download.path,"hgnc_complete_set.txt.gz",sep="/")))
           download.file(hgnc,paste(download.path,"hgnc_complete_set.txt.gz",sep="/"),method="auto")
       hgnc <- paste(download.path,"hgnc_complete_set.txt.gz",sep="/")
    }
    if(substr(hgnc,nchar(hgnc)-1,nchar(hgnc)) == "gz"){
        hgnc <- read.delim(gzfile(hgnc))
        }else{
            hgnc <- read.delim(hgnc)
    }       
   # HPO
       if(is.null(HPO.disease)){
          HPO.disease.check <- pheno_extract_HPO(keyword= keyword)
          HPO.disease <- as.character(unique(HPO.disease.check[grep("OMIM",HPO.disease.check[,1]),1]))
       }

   #start to search
       if(!is.null(keyword)){
          morbidmap.d <- morbidmap[grep_split(keyword,morbidmap[,1]),]
          }else if((is.null(keyword)& !is.null(HPO.disease)) | (is.null(keyword) & !is.null(genelist))){
            morbidmap.d <- c()  
            }else{
              morbidmap.d <- morbidmap
       }           

       if(!is.null(HPO.disease)){
          for(i in HPO.disease){
              i.trim <- unlist(strsplit(i,"OMIM:"))[2]
              morbidmap.d <- rbind(morbidmap.d, morbidmap[grep_split(keyword=i.trim,morbidmap[,1]),])
          }
       }
#add the genelist,may capture other phenotypes
      genes.omim <- lapply(as.character(morbidmap[,2]),function(x) unlist(strsplit(gsub(" ", "",x),",")))
      for(j in as.character(genelist)){
         m.j <- morbidmap[unlist(lapply(genes.omim,function(x) is.element(j,x))),]
         morbidmap.d <- rbind(morbidmap.d, m.j)
      }
       morbidmap.d <- unique(morbidmap.d)
       

       gene.3 <- grep(" (3)",morbidmap.d[,1],fixed = TRUE)  ;# 3 represent the molecular basis for the disorder is known; 
                                                       # a mutation has been found in the gene.
       genes <- unique(morbidmap.d[gene.3,"gene.mim.no"])
       hgnc.extract <- hgnc[is.element(hgnc[,34],morbidmap.d[,3]),]
       missing <- as.character(setdiff(morbidmap.d[,3],hgnc.extract[,34]))
       morbidmap.d$Approved.Symbol <- morbidmap.d$Approved.Name <- "missing"
       for(i in hgnc.extract[,34]){
          morbidmap.d[grep(i,morbidmap.d[,3]),c("Approved.Symbol")] <- as.character(hgnc.extract[grep(i,hgnc.extract[,34]),c("Approved.Symbol")])
          morbidmap.d[grep(i,morbidmap.d[,3]),c("Approved.Name")] <- as.character(hgnc.extract[grep(i,hgnc.extract[,34]),c("Approved.Name")])
       }

#check the missing gene symbols
       missing.mimnos <- as.character(unique(morbidmap.d[morbidmap.d$Approved.Symbol == "missing",3]))
       for(g in missing.mimnos){
         hgnc.g <- as.character(hgnc[grep(g,hgnc[,34]),"Approved.Symbol"])
         if(length(hgnc.g) ==1 ){
            missing <- setdiff(missing,g)
            morbidmap.d[grep(g,morbidmap.d[,3]),"Approved.Symbol"] <- as.character(hgnc[grep(g,hgnc[,34]),"Approved.Symbol"])
            morbidmap.d[grep(g,morbidmap.d[,3]),c("Approved.Name")] <- as.character(hgnc[grep(g,as.character(hgnc[,34])),"Approved.Name"])
         }
       }
       if(length(missing)>0){
           print(paste(paste(missing,collapse=","),ifelse(length(missing) ==1," does"," do")," NOT contained in HGNC database.\n"))
       }


#check.in(): to check whether a keyword is in a string 
    check.in <- function(y){ 
         ifelse(length(grep(toupper(keyword),as.character(y)))>0,"yes","no")
    }

   if(type != "gene"){
      variants <- c()
      for(j in genes){
          print(j)
          urls <- paste("http://api.europe.omim.org/api/entry?mimNumber=",j,"&include=all&apiKey=",omim.apiKey,sep="")
          doc <- xmlTreeParse(urls,useInternalNodes = TRUE)       
          nodes <- getNodeSet(doc, "//entry")
          lists <- nodesToList(nodes)
      
      #gene description
          mimNumber.j  <- j     
          geneName.j <- lists[[1]]$geneMap$geneName
          chromosome.j <- lists[[1]]$geneMap$chromosome
          cytoLocation.j <- lists[[1]]$geneMap$cytoLocation
          Approved.Symbol.j <- unique(morbidmap.d[morbidmap.d[,3]==j,6])
     
      #phenotypeMapList
          phenotypeMapList.j <- lists[[1]]$geneMap$phenotypeMapList
          phenotypeMap.j <- names.inheri <- c()
          if(length(phenotypeMapList.j) > 0 ){
              for(m in 1:length(phenotypeMapList.j)){
                 phenotype.m <- phenotypeMapList.j[[m]]$phenotype
                 phenotypeInheritance.m <- phenotypeMapList.j[[m]]$phenotypeInheritance
                 phenotypeInheritance.m <- ifelse(is.null(phenotypeInheritance.m),"",phenotypeInheritance.m)
                 phenotypeMap.j <- c(phenotypeMap.j,phenotypeInheritance.m)
                 names.inheri <- c(names.inheri, phenotype.m)
              }
              names(phenotypeMap.j) <- toupper(names.inheri)
          }
          
      #variants
          variants.exists.j <- lists[[1]]$allelicVariantExists
        if(!is.null(variants.exists.j)){
          if(variants.exists.j == "true") {
              allelic.variant.j <- lists[[1]]$allelicVariantList
              variants.j <- matrix(,ncol=9,nrow=length(allelic.variant.j))
              colnames(variants.j) <- c("number","variants.ID","status","phenotype","mutations","text","dbSNPs","clinvarAccessions","Inheritance")
              for(i in 1:length(allelic.variant.j)){
                 variants.j.i <- allelic.variant.j[[i]]
                 variants.j[i,"number"] <- variants.j.i$number                 
                 variants.j[i,"status"] <- variants.j.i$status
                 variants.j[i,"phenotype"] <- variants.j.i$name
                 variants.j[i,"Inheritance"] <- phenotypeMap.j[variants.j.i$name]
                 if( variants.j.i$status != "moved" & variants.j.i$status != "removed"){
                    mut.j.i <- variants.j.i$mutations
                    variants.j[i,"mutations"] <- ifelse(length(mut.j.i)==1, mut.j.i, ifelse(length(mut.j.i)>1,paste(mut.j.i,collapse=";"),""))
                    text.j.i <- variants.j.i$text
                    if(length(grep("\n",text.j.i))>0) {
                       text.j.i <- unlist(strsplit(text.j.i,"\n"))
                       text.j.i <- text.j.i[text.j.i != ""]
                       text.j.i <- paste(text.j.i,collapse="..")
                    }
                    if(is.null(text.j.i))  text.j.i <- ""
                    variants.j[i,"text"] <- text.j.i
                    variants.j[i,"dbSNPs"] <- ifelse(!is.null(variants.j.i$dbSnps),variants.j.i$dbSnps,"")
                    variants.j[i,"clinvarAccessions"] <- ifelse(!is.null(variants.j.i$clinvarAccessions),variants.j.i$clinvarAccessions,"")
                 }
              }
              if(i >1){
                 variants <- rbind(variants,cbind(mimNumber.j,Approved.Symbol.j,geneName.j,chromosome.j,cytoLocation.j,variants.j))
                 }else if(i == 1) {
                   variants <- rbind(variants,c(mimNumber.j,Approved.Symbol.j,geneName.j,chromosome.j,cytoLocation.j,variants.j))
              }   
              }}else{
                   variants <- rbind(variants,c(mimNumber.j,Approved.Symbol.j,geneName.j,chromosome.j,cytoLocation.j,rep("",9)))
          }       
      }    

      phenotype.check.all <- unlist(lapply(variants[,"phenotype"], check.in))
      variants <- cbind(variants, phenotype.check.all)
      colnames(variants) <- c("mimNumber.gene","Approved.Symbol","GeneName","Chromosome","cytoLocation","Variants.number","variants.ID","status","Phenotype",
                             "mutations","text", "dbSNPs","clinvarAccessions", "Inheritance","phenotype.check")
      variants[,"variants.ID"] <- paste("OMIM Allelic Variant:", variants[,"mimNumber.gene"], ".", variants[,"Variants.number"],sep="")
   }      
   if(type == "both"){
      extract <- list(morbidmap.d,variants)
      names(extract) <- c("morbidmap","variants")
      }else if (type == "gene"){
           extract <- morbidmap.d
           }else if(type == "variant"){
            extract <- variants
   }         
   return(extract)
}
