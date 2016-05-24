extract_pubmed <-
function(query, keyword,localPDB = paste(getwd(),"localPDB",sep="/")){ 
    ## extract the RESULTS from abstract
    abs_trim <- function(x){
       # x = pubmed_abs[1]
        if(length(grep("METHODS:",x)) > 0) {
           x.trim <- unlist(strsplit(x,"METHODS:"))[2]
           }else if(length(grep("METHOD:",x)) > 0) {
               x.trim <- unlist(strsplit(x,"METHOD:"))[2]
               }else {
                 x.trim <- x  
        }
        if(length(grep("CONCLUSION",x.trim)) > 0) 
              x.trim <- unlist(strsplit(x.trim,"CONCLUSION*"))[1]     
        return(x.trim)  
    }
    
    # extract the conclusion from abstract
    abs_conclusion <- function(x){
       # x = pubmed_abs[1]
        if(length(grep("CONCLUSIONS:",x)) > 0) {
           x.trim <- unlist(strsplit(x,"CONCLUSIONS:"))[2]
           }else if(length(grep("CONCLUSION:",x)) > 0) {
                 x.trim <- unlist(strsplit(x,"CONCLUSION:"))[2]
                 }else{
                 x.trim <- x  
        }
        return(x.trim)  
    }
    
# transform the gene alias to gene symbol
gs2hgnc = function(gene){
   #gene = "A1BGAS"
   gene = toupper(gene)
   gene.approved = gene
   app = hgnc[hgnc$Approved.Symbol==gene,]
   pre.symbol = hgnc[grep(gene,hgnc$Previous.Symbols),]
   synonyms = hgnc[grep(gene,hgnc$Synonyms),]
   if(nrow(app)==1){
      gene.approved = as.character(app$Approved.Symbol)
      }else if(nrow(pre.symbol)==1){
         pre.symbol.trim = as.character(pre.symbol$Previous.Symbols)
         pre.symbol.trim = unlist(strsplit(pre.symbol.trim,", "))
         if(length(pre.symbol.trim[pre.symbol.trim==gene]) == 1){
            gene.approved = as.character(pre.symbol$Approved.Symbol)
            }
          }else if(nrow(pre.symbol) > 1){
            pre.symbol.trim = as.character(pre.symbol$Previous.Symbols)
            for(j in 1:length(pre.symbol.trim)){
                  pre.symbol.trim.j = unlist(strsplit(pre.symbol.trim[j],", "))
                  if(length(intersect(pre.symbol.trim.j,gene)) > 0){
                  if(intersect(pre.symbol.trim.j,gene) == gene){
                    pre.symbol.j = pre.symbol[j,]
                    gene.approved = as.character(pre.symbol.j$Approved.Symbol)
                  }}
               }
            }else if(nrow(synonyms) >= 1){
                synonyms.trim = as.character(synonyms$Synonyms)
                for(j in 1:length(synonyms.trim)){
                  synonyms.trim.j = unlist(strsplit(synonyms.trim[j],", "))
                  if(length(intersect(synonyms.trim.j,gene)) > 0){
                  if(intersect(synonyms.trim.j,gene) == gene){
                   synonyms.j = synonyms[j,]
                    gene.approved = as.character(synonyms.j$Approved.Symbol)
                  }}
               }               
            }
   if(nrow(synonyms) == 0 & nrow(app)==0 & nrow(pre.symbol)==0) { gene.approved = paste("missing(",gene,")",sep="") }        
   return(gene.approved)
}

    # x is a string, title, background, conclusion, results etc
    # pheno_capture_abs(keyword,pubmed_title[91])
    pheno_capture_abs <- function(keyword,x){
         # x <- pubmed_title[16]
         x.trim <- unlist(strsplit(x," with |\\.|\\:| cause. | associated | mapped 
                    |\\,| due | to | by | was | were | that | of | in | on | is 
                    | are | the | a | an | for | identified | by | using | .ing 
                    | susceptibility| had | featuring | and therefore"))
         x.trim <- unlist(strsplit(x.trim,"with |\\.|caus. |\\,|to |by |was |were 
               |that |of | in | is| are|for |the | susceptibility |had |causing"))
         pheno.extract <-  paste(unique(x.trim[c(grep_split(keyword,x.trim),
                                grep("congenital",x.trim),
                                grep("deficiency",x.trim),
                                grep("susceptibility",x.trim),
                                grep("X-linked",x.trim),
                                grep("dominant",x.trim),
                                grep("autosomal",x.trim),
                                grep("dominant",x.trim),
                                grep("recessive",x.trim),
                                grep("dysplasia",x.trim),
                                grep("familial",x.trim),
                                grep("dystrophy",x.trim),
                                grep("disorder",x.trim),
                                grep("ataxia",x.trim),
                                grep("hereditary",x.trim),
                                grep("Mitochondrial",x.trim),
                                grep("hypoplasia",x.trim),
                                grep("Mitochondrial",x.trim),
                                grep("syndrome",x.trim))]),collapse=";")
         if(pheno.extract == ""){
            x.trim <- unlist(strsplit(x," with "))[2]
            pheno.extract <- unlist(strsplit(x.trim,"\\.| by "))[1]
            pheno.extract <- ifelse(length(grep("mutation",pheno.extract)) > 0, 
                                        pheno.extract[-1],pheno.extract)
            if(is.na(pheno.extract)){
               if(length(grep("X.linked",x)) > 0 ){
                  x.trim <- unlist(strsplit(x," X.linked "))[2]
                  pheno.extract <- paste("X-linked", 
                          unlist(strsplit(x.trim,"\\.| by | in | with "))[1],sep=" ")
               }
               if(length(grep("Y.linked",x)) > 0 ){
                  x.trim <- unlist(strsplit(x,"Y.linked "))[2]
                  pheno.extract <- paste("Y-linked", 
                       unlist(strsplit(x.trim,"\\.| by | in | with "))[1],sep=" ")
               }
               if(length(grep("autosomal dominant",x)) > 0 ){
                  x.trim <- unlist(strsplit(x,"autosomal dominant"))[2]
                  pheno.extract <- paste("autosomal dominant",
                      unlist(strsplit(x.trim,"\\.| by | in | with "))[1],sep=" ")
               }
               if(length(grep("autosomal recessive",x)) > 0 ){
                  x.trim <- unlist(strsplit(x,"autosomal recessive"))[2]
                  pheno.extract <- paste("autosomal recessive",
                      unlist(strsplit(x.trim,"\\.| by | in | with "))[1],sep=" ")
               }
            }
         }   
         pheno.extract <- stri_replace_all_regex(pheno.extract,
              '\\:$|\\]|\\[|^ and|^ | $', '')
         if(length(pheno.extract) >1)
            pheno.extract <-  paste(pheno.extract,collapse=";")
          return(pheno.extract)
    }
    
    ## capture the gene and mutation from an abstract
    gene_var_abs <- function(x){
        # x = abs_trim(pubmed_abs[93])
        x <- stri_replace_all_regex(x, " &gt; ",  "&gt;") 
        x <- stri_replace_all_regex(x, "-&gt;",  "&gt;") 
        x.trim <- unlist(strsplit(x,"( )|\\)|\\(|\\:|\\,|\\[|\\]|\\-associated"))
        x.trim <- stri_replace_all_regex(x.trim,"\\;$","")
         n_c. <- grep("(^c\\.[0-9]+$|^c\\.$)",x.trim)
         if(length(n_c.) > 0){
           for(i in 1:length(n_c.)){
             n_c.i <- n_c.[i]
             x.trim <- c(x.trim,paste(x.trim[n_c.i], x.trim[n_c.i + 1],sep=""))
             x.trim <- x.trim[-(n_c.i:(n_c.i+1))]
           }  
         }
         n_p. <- grep("^p\\.$",x.trim)
         if(length(n_p.) > 0){
           for(i in 1:length(n_p.)){
              n_p.i <- n_p.[i]
              x.trim <- c(x.trim,paste("p.", x.trim[n_p.i + 1],sep=""))
           }   
              x.trim <- x.trim[-c(n_p.,n_p.+1)]
         }
         n_p.space <- grep("^p\\.[A-Z][0-9]+$",x.trim)
         if(length(n_p.space) > 0){
           for(i in 1:length(n_p.space)){
              n_p.space.i <- n_p.space[i]
              x.trim <- c(x.trim,paste(x.trim[n_p.space.i], x.trim[n_p.space.i + 1],sep=""))
           }   
              x.trim <- x.trim[-c(n_p.space,n_p.space+1)]
         }
             
         x.trim <- stri_replace_all_regex(x.trim, '\\.$', '')
         muts.dna.1 <- x.trim[grep("^c\\..",x.trim)]
         muts.dna.2 <- x.trim[grep("[ATCG]\\&gt\\;[ATCG]",x.trim)]   
         muts.dna.3 <- x.trim[grep("^IVS.",x.trim)]  
         muts.dna.4 <- x.trim[grep("[0-9]+ins[ATCG]",x.trim)]  
         muts.dna.5 <- x.trim[grep("rs[0-9]+",x.trim)] 
         muts.dna.6 <- x.trim[grep("[0-9]+del[ATCG]",x.trim)]                   
         muts.dna.7 <- x.trim[grep("[0-9]+\\_[0-9]dup*.",x.trim)]                   
         muts.dna <- unique(c(muts.dna.1,muts.dna.2,muts.dna.3, muts.dna.4,muts.dna.5,muts.dna.6,muts.dna.7))
         
         #search the AA mutation
         ## p.* format in the first
         muts.pt <- unique(x.trim[grep("^p\\.",x.trim)])
         muts.pt.trim <- c()
         if(length(muts.pt) >= 1){
           for(i in 1:length(muts.pt)){
              muts.pt.1 <- muts.pt[i]
              muts.pt.trim.i = c()
              if(length(grep("^p\\.[A-Z][0-9]+[A-Z][0-9]+",muts.pt)) ==1){              
                  muts.pt.1.split <- unlist(strsplit(muts.pt.1,"p\\.|[0-9]+|\\_|[0-9]+|dup|del"))
                  muts.pt.hgvs.1 <- muts.pt.1.split[is.element(muts.pt.1.split,aa.table[,4])]
                  muts.pt.trim.i <- stri_replace_all_regex(muts.pt.1,paste("^p\\.",muts.pt.hgvs.1[1],sep=""),paste("p.",as.character(unique(aa.table[aa.table[,4] == muts.pt.hgvs.1[1],3])),sep=""))
                  muts.pt.trim.i <- stri_replace_all_regex(muts.pt.trim.i,paste("_",muts.pt.hgvs.1[2],sep=""),paste("_",as.character(unique(aa.table[aa.table[,4] == muts.pt.hgvs.1[2],3])),sep=""))              
                  muts.pt.trim <- c(muts.pt.trim,muts.pt.trim.i)
              }else if(length(grep("del",muts.pt.1,ignore.case=T))==1){
                  muts.pt.1.split <- unlist(strsplit(muts.pt.1,"del"))
                  if(is.element(muts.pt.1.split[2],aa.table[,4])){
                      muts.pt.1 <- stri_replace_all_regex(muts.pt.1,muts.pt.1.split[2],"")
                      muts.pt.trim.i <- stri_replace_all_regex(muts.pt.1,"^p\\.",paste("p.",as.character(unique(aa.table[aa.table[,4] == muts.pt.1.split[2],3])),sep=""))
                      muts.pt.trim <- c(muts.pt.trim,muts.pt.trim.i)
                     }else{muts.pt.trim <- ""}                      
                  }else{
                       if(length(grep("fs",muts.pt.1,ignore.case=T))==1)
                          muts.pt.1 <- unlist(strsplit(muts.pt.1,"fsX|fs|fs\\*"))[1]
                       muts.pt.hgvs <- unlist(strsplit(muts.pt.1,"p\\.|[0-9]+"))
                       muts.pt.hgvs <- muts.pt.hgvs[muts.pt.hgvs != ""]
                       muts.pt.hgvs.1 <- muts.pt.hgvs[is.element(muts.pt.hgvs,aa.table[,4])]
                       if(length(muts.pt.hgvs.1) == 1){
                          if(length(grep("\\*$",muts.pt.1,ignore.case=T))==1){ 
                             muts.pt.hgvs.1_3_1 <- as.character(unique(aa.table[aa.table[,4] == muts.pt.hgvs.1[1],3]))
                             muts.pt.trim.i <- stri_replace_all_regex(muts.pt.1,paste("p.",muts.pt.hgvs.1[1],sep=""),paste("p.",muts.pt.hgvs.1_3_1,sep=""))
                             } else{
                                muts.pt.trim.i = c()
                          }                 
                       }
                       
                       if(length(muts.pt.hgvs.1) == 2){
                          muts.pt.hgvs.1_3_1 <- as.character(unique(aa.table[aa.table[,4] == muts.pt.hgvs.1[1],3]))
                          muts.pt.hgvs.1_3_2 <- as.character(unique(aa.table[aa.table[,4] == muts.pt.hgvs.1[2],3]))
                          muts.pt.trim.i <- stri_replace_all_regex(muts.pt.1,paste("p.",muts.pt.hgvs.1[1],sep=""),paste("p.",muts.pt.hgvs.1_3_1,sep=""))
                          if(length(grep("fs",muts.pt[1],ignore.case=T))==1){
                            muts.pt.trim.i <- stri_replace_all_regex(muts.pt.trim.i,muts.pt.hgvs.1[2],paste(muts.pt.hgvs.1_3_2,'fs',sep=""))
                            }else{
                              muts.pt.trim.i <- stri_replace_all_regex(muts.pt.trim.i,paste(muts.pt.hgvs.1[2],'$',sep=""),muts.pt.hgvs.1_3_2)                            
                          }                     
                       }
                     }  
               ## check the 3-letters AA
              muts.pt.hgvs <- unlist(strsplit(muts.pt.1,"p\\.|[0-9]+"))
              muts.pt.hgvs <- muts.pt.hgvs[muts.pt.hgvs != ""]
              if(is.element(muts.pt.hgvs[1],aa.table[,3])) muts.pt.trim.i <- muts.pt[i]
              if(length(muts.pt.trim.i) > 0){
                         muts.pt.trim <- c(muts.pt.trim,muts.pt.trim.i) 
                         rm(muts.pt.trim.i)    
              }                                                   
           } 
         }

         
         ## if can not find, then AA+numeric+AA
  ##       if(length(muts.pt) == 0 ){
            muts.pt.2 <- unique(x.trim[grep("^[A-Z]*.*[0-9]+[A-Z]",x.trim)])
          #  muts.pt.2 <- muts.pt.2[-grep("^\\p.",muts.pt.2)]
            muts.pt.2 <- setdiff(muts.pt.2, muts.pt.2[grep("^D[0-9]+S[0-9]+",muts.pt.2)])
            muts.pt.2 <- setdiff(muts.pt.2, muts.pt.2[grep("^D[XY]S[0-9]+",muts.pt.2)])
            muts.pt <- union(muts.pt,muts.pt.2)
            ## check the AA mutation
            if(length(muts.pt.2 ) > 0){
             for( i in muts.pt.2){
                aa.check = unlist(strsplit(i,'[0-9]+'))
                aa.check = aa.check[aa.check != ""]
                aa.check = unique(is.element(aa.check,aa_ab))
                aa.check = aa.check[order(aa.check)]
                if(!aa.check[1])
                   muts.pt = setdiff(muts.pt,i)
             }}
         if(length(muts.pt.2) >= 1){
           for(i in 1:length(muts.pt.2)){
              muts.pt.1 <- muts.pt.2[i]
              muts.pt.trim.i = c()
                       muts.pt.hgvs <- unlist(strsplit(muts.pt.1,"[0-9]+"))
                       muts.pt.hgvs.1 <- muts.pt.hgvs[is.element(muts.pt.hgvs,aa.table[,4])]
                       muts.pt.trim.i = c()
                       if(length(muts.pt.hgvs.1) == 2){
                          muts.pt.hgvs.1_3_1 <- as.character(unique(aa.table[aa.table[,4] == muts.pt.hgvs.1[1],3]))
                          muts.pt.hgvs.1_3_2 <- as.character(unique(aa.table[aa.table[,4] == muts.pt.hgvs.1[2],3]))
                          muts.pt.trim.i <- stri_replace_all_regex(muts.pt.1,muts.pt.hgvs.1[1],paste("p.",muts.pt.hgvs.1_3_1,sep=""))
                          muts.pt.trim.i <- stri_replace_all_regex(muts.pt.trim.i,paste(muts.pt.hgvs.1[2],'$',sep=""),muts.pt.hgvs.1_3_2)
                          }
                      if(length(muts.pt.trim.i) > 0){
                         muts.pt.trim <- c(muts.pt.trim,muts.pt.trim.i) 
                         rm(muts.pt.trim.i)    
                      }               
           }  
         }    
         ## check whether the mutation is AA+nmeric+*, premature termanition
            muts.pt.3 <- unique(x.trim[grep("^[A-Z]*.*[0-9]+\\*$",x.trim)])
            if(length(muts.pt.3 ) > 0){
              for(i in muts.pt.3 ){
                if(length(grep(i,muts.pt)) == 0)
                  muts.pt <- union(muts.pt,muts.pt.3)
              }    
            }
            if(length(muts.pt.3) >0){
              for(i in 1:length(length(muts.pt.3))){
                muts.pt.1 <- muts.pt.3[i]
                muts.pt.trim.i = c()
                muts.pt.hgvs <- unlist(strsplit(muts.pt.1,"[0-9]+"))
                muts.pt.hgvs.1 <- muts.pt.hgvs[is.element(muts.pt.hgvs,aa.table[,4])]
                if(length(grep("\\*$",muts.pt.1,ignore.case=T))==1){ 
                   muts.pt.hgvs.1_3_1 <- as.character(unique(aa.table[aa.table[,4] == muts.pt.hgvs.1[1],3]))
                   muts.pt.trim.i <- stri_replace_all_regex(muts.pt.1,muts.pt.hgvs.1[1],paste("p.",muts.pt.hgvs.1_3_1,sep=""))
                   } else{
                      muts.pt.trim.i = c()
                }
                muts.pt.trim = c(muts.pt.trim,muts.pt.trim.i) 
                rm(muts.pt.trim.i)
              }                 
            }
            
            
         ## remove the STRs, such as D1S498
            if(length(muts.pt) > 0){
                muts.pt <- setdiff(muts.pt, muts.pt[grep("^D[0-9]+S[0-9]+",muts.pt)])
                muts.pt <- setdiff(muts.pt, muts.pt[grep("^D[XY]S[0-9]+",muts.pt)])
            }
     ##    }
         ## last, search AA full name mutation
         if(length(muts.pt) == 0 ){
            muts.no <- unlist(lapply(aa_full,function(x.aa) agrep(x.aa,x.trim,value=F)))
            muts.no <- muts.no[order(muts.no)]
            if(length(muts.no) > 0 ){
               muts.pt.4 <- x.trim[muts.no]
               muts.pt <- muts.pt.4
               for(i in 1:length(muts.pt.4)){
                  muts.pt.4.trim <- as.character(aa.table[aa.table[,2]==muts.pt.4[i],3][1])
                  muts.pt.trim <-  c(muts.pt.trim, muts.pt.4.trim)
               } 
            }            
         }
         muts.pt <- unique(muts.pt)
    ## search genes
         genes <- x.trim[grepl("^[A-Z]*.*[A-Z0-9]", x.trim)]
         genes.symbol <- unique(genes[is.element(genes,Approved.Symbol)])
         genes.symbol <- genes.symbol[!is.element(genes.symbol,c("A","C","G","T","I","III","II","IV",
                   "V","VI","VII","VIII","IVV","VV","DNA","RNA","HDR","LOD","WT","MRI"))]
         genes.alias <- unique(genes[is.element(genes,gene.Alias)])
         genes.alias <- genes.alias[nchar(genes.alias)>2]
         genes.alias <- genes.alias[!is.element(genes.alias,c("A","C","G","T","I","III","II","IV",
                   "V","VI","VII","VIII","IVV","VV","DNA","RNA","HDR","LOD","WT","BLAST"))]
         if(length(genes.symbol) == 0)  
              genes.symbol = ""
         if(length(genes.alias) == 0) { 
              genes.alias = ""
              }else{
                   genes.s <- unlist(lapply(genes.alias,gs2hgnc))
                   genes.symbol <- union(genes.symbol,genes.s)
         }  
         
         if(length(muts.dna) == 0) {
            muts.dna.trim <- ""
            }else{
              muts.dna <- stri_replace_all_regex(muts.dna, '\\,$', '')
              muts.dna.trim <- stri_replace_all_regex(muts.dna, '\\&gt\\;', '>')  
         }
         if(length(genes.symbol)>1) 
                   genes.symbol <- paste(unique(genes.symbol),collapse=", ")  
         if(length(genes.alias)>1) 
                   genes.alias <- paste(unique(genes.alias),collapse=", ")  
         if(length(muts.dna.trim)>1) 
                   muts.dna.trim <- paste(unique(muts.dna.trim),collapse=", ")  
         if(length(muts.pt) == 0 ) {
            muts.pt <- ""
            }else if(length(muts.pt)>1) {
               muts.pt <- paste(unique(muts.pt),collapse=", ") 
         }       
         if(length(muts.pt.trim) == 0 ) {
            muts.pt.trim <- ""
            }else if(length(muts.pt.trim)>1) {
               muts.pt.trim <- paste(unique(muts.pt.trim),collapse=", ") 
         }       
         caps = c(genes.symbol,genes.alias,muts.dna.trim,muts.pt,muts.pt.trim)    
         return(caps)
    }

## input HGNC dataset
    if(file.exists(localPDB)){
         if(file.exists(paste(localPDB,"hgnc_complete_set.txt.gz",sep="/"))){
             hgnc <- paste(localPDB,"hgnc_complete_set.txt.gz",sep="/")
             }else{
                 hgnc <- NULL
         }        
        }else{
             hgnc <- NULL
    }     

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

  ##  hgnc <- read.delim(gzfile("/public/home/czf/project/rare.disease/localPDB/hgnc_complete_set.txt.gz"))
    Approved.Symbol <- as.character(hgnc$Approved.Symbol)
    gene.Synonyms <- unlist(lapply(as.character(hgnc$Synonyms),function(x) str_trim(unlist(strsplit(x,",")))))
    gene.Previous.Symbols <- unlist(lapply(as.character(hgnc$Previous.Symbols),function(x) str_trim(unlist(strsplit(x,",")))))
    gene.Alias <- c(gene.Synonyms,gene.Previous.Symbols)
#    aa <- NULL
#    data(aa,package="VarfromPDB")
    aa.table <- aa
    aa_full = unique(as.character(aa.table[,2]))
    aa_ab = as.character(unique(as.matrix(aa.table[,3:4])))
    
    ## search in PubMed
    pubmed_search <- EUtilsSummary(query, type="esearch",db = "pubmed",retmax=30000)
    pubmed_records <- EUtilsGet(pubmed_search)
    years <- YearPubmed(pubmed_records)
    authors <- Author(pubmed_records)
    first_author <- unlist(lapply(authors,function(x) paste(as.character(x[1,1:2]),collapse=" ")))
    country <- Country(pubmed_records)
    pubmed_abs <- AbstractText(pubmed_records)
    pubmed_Affiliation <- Affiliation(pubmed_records)
    pubmed_title <- ArticleTitle(pubmed_records)
    pubmed_PMID <- PMID(pubmed_records)
    pubmed_abs.trim <- unlist(lapply(pubmed_abs,abs_trim))
    pubmed_conclusion <- unlist(lapply(pubmed_abs,abs_conclusion))
    pubmed_captures.1 <- t(matrix(unlist(lapply(pubmed_title,gene_var_abs)),nrow=5))
    pubmed_captures <- t(matrix(unlist(lapply(pubmed_abs,gene_var_abs)),nrow=5))
    pubmed_captures[pubmed_captures.1[,1] != "",1] =  pubmed_captures.1[pubmed_captures.1[,1] != "",1]
#    pubmed_captures <- t(matrix(unlist(lapply(pubmed_abs.trim,gene_var_abs)),nrow=5))
    phenotype_pubmed <- unlist(lapply(pubmed_title, function(x) pheno_capture_abs(keyword,x)))
    phenotype_pubmed[is.na(phenotype_pubmed)] <- unlist(lapply(pubmed_conclusion[is.na(phenotype_pubmed)], function(x) pheno_capture_abs(keyword,x)))
    pubmed_captures <- cbind(phenotype_pubmed,pubmed_captures,pubmed_title,years,first_author,country,pubmed_PMID)
    colnames(pubmed_captures) <- c("Phenotype","Approved.Symbol", "Genes.possible","c.DNA_change", "p.change","p.change.HGVS", "Article_Title", "Year", "First_author", "Country", "PMID")
    print(c(length(phenotype_pubmed),length(pubmed_title),length(years),length(first_author),length(country),length(pubmed_PMID)))
    print(dim(pubmed_captures))
    return(list(pubmed_captures,pubmed_abs))    
}
