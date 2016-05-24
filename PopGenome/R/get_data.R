get_data <- function(matr,include.unknown=FALSE,gff=FALSE,FAST,SNP.DATA){

reverse.codons <- FALSE

# poly <- get.polymorph(matr)
# matr <- matr[,poly,drop=FALSE]
# cat("Suche nach polymorphen Stellen ist fertig !")
# GLOBAL is an evironment

if(FAST){
 
  
  if(include.unknown){
  bialpos <- .Call("polyCinclude",matr)
  }else{  
  bialpos <- .Call("polyC",matr)
  }


#bialpos <- as.logical(bialpos)
#RETURNLISTE <- list(n_site=n_site,transitions=matrix_sv, biallelic.matrix=matrix_pol,
#biallelic.sites=matrix_pos,matrix_codonpos=matrix_codonpos,n.singletons=unic,totalmis=NaN,s_sites=NaN,mvbp=mvbp
#,trans.transv.ratio=(transitions/transversions),n.valid.sites=algsites,n.biallelic.sites=bial_sites,
#polyallelic.sites=mhitbp,n.nucleotides=sum_sam,biallelic.compositions=TCGA,ROUGH=erg,matrix_freq=matrix_freq,syn=syn,
#nonsyn=nonsyn,synonymous=!as.logical(synnonsyn),biallelic.substitutions=subst,minor.alleles=mutations,codons=Codons,
#CodingSNPS=CodingSNPS,UTRSNPS=UTRSNPS,IntronSNPS=IntronSNPS,ExonSNPS=ExonSNPS,GeneSNPS=GeneSNPS,Coding_region_length=Coding_region_length,
#UTR_region_length=UTR_region_length,Intron_region_length=Intron_region_length,Exon_region_length=Exon_region_length, #Gene_region_length=Gene_region_length, sites.with.gaps=gaps,sites.with.unknowns=unknown) 

 
biallelic.sites    <- which(bialpos==1)
polyallelic.sites  <- which(bialpos==4)

#FIXME
gaps            <- which(bialpos==2)
unknowns        <- which(bialpos==3)


if(length(biallelic.sites)==0){
# print("No biallelic positions !")
return(NA)
}


 SUBMAT          <- matr[,biallelic.sites,drop=FALSE]

 ## very fast ---------------
 if(include.unknown){
 res          <- .Call("makeBialMatrixinclude",SUBMAT)
 Bial.Mat     <- res[[1]]
 Bial.Mat[Bial.Mat==-1] <- NaN
 }else{
 res           <- .Call("makeBialMatrix",SUBMAT)
 Bial.Mat      <- res[[1]]
 }



 transitions             <- as.logical(res[[2]])
 biallelic.substitutions <- res[[3]]
 rownames(biallelic.substitutions) <- c("minor","major")
 n.transitions           <- sum(transitions)
 n.transversions         <- length(biallelic.sites) - n.transitions
 tt.ratio                <- n.transitions/n.transversions 



## GFF-File
# GFF
#if(is.list(gff)){
#features <- parse_gff(gff)
#}

## GFF
if(is.list(gff)){

features <- gff

   if(length(features$Gene)>0){
    Gene_region        <- as.vector(unlist(apply(features$Gene,1,function(x){return(x[1]:x[2])})))
    GeneSNPS           <- is.element(biallelic.sites,Gene_region)
    Gene_region_length <- length(Gene_region)
    #IntronSNPS    <- matrix_pos[IntronSNPS]
   }else{GeneSNPS<-NaN;Gene_region_length <- 0} 

   if(length(features$Intron)>0){
    Intron_region <- as.vector(unlist(apply(features$Intron,1,function(x){return(x[1]:x[2])})))
    IntronSNPS    <- is.element(biallelic.sites,Intron_region)
    Intron_region_length <- length(Intron_region)
    #IntronSNPS    <- matrix_pos[IntronSNPS]
   }else{IntronSNPS<-NaN;Intron_region_length <- 0}

   if(length(features$UTR)>0){	
    UTR_region        <- as.vector(unlist(apply(features$UTR,1,function(x){return(x[1]:x[2])})))
    UTRSNPS           <- is.element(biallelic.sites,UTR_region)
    UTR_region_length <- length(UTR_region) 
    #UTRSNPS       <- matrix_pos[UTRSNPS]
   }else{UTRSNPS<-NaN;UTR_region_length <- 0}

   if(length(features$Exon)>0){
    Exon_region        <- as.vector(unlist(apply(features$Exon,1,function(x){return(x[1]:x[2])})))
    ExonSNPS           <- is.element(biallelic.sites,Exon_region)
    Exon_region_length <- length(Exon_region)
   }else{ExonSNPS<-NaN;Exon_region_length<-0}
   
   if(length(features$Coding)>0){ 
    Coding_region  <- as.vector(unlist(apply(features$Coding,1,function(x){return(x[1]:x[2])})))
    CodingSNPS     <- is.element(biallelic.sites,Coding_region)
    Coding_region_length <- length(Coding_region)
   # CodingSNPS     <- biallelic.sites[CodingSNPS]
   # size           <- length(CodingSNPS)
   }else{CodingSNPS<-NaN;Coding_region_length <- 0}

 }else{CodingSNPS<- NaN;IntronSNPS <- NaN; UTRSNPS <- NaN;ExonSNPS <- NaN;GeneSNPS<-NaN;
Coding_region_length<-NaN;Intron_region_length<-NaN;UTR_region_length<-NaN;Exon_region_length<-NaN;Gene_region_length<-NaN}
  

return(list(biallelic.matrix=Bial.Mat,biallelic.sites=biallelic.sites,polyallelic.sites=polyallelic.sites,
sites.with.gaps=gaps,sites.with.unknowns=unknowns,
transitions=transitions,n.valid.sites=NaN,synonymous=rep(NaN,length(biallelic.sites)),
trans.transv.ratio=tt.ratio,biallelic.substitutions=biallelic.substitutions,CodingSNPS=CodingSNPS,UTRSNPS=UTRSNPS,IntronSNPS=IntronSNPS,
ExonSNPS=ExonSNPS,GeneSNPS=GeneSNPS,Coding_region_length=Coding_region_length,
Gene_region_length=Gene_region_length, Exon_region_length=Exon_region_length, Intron_region_length=Intron_region_length,UTR_region_length=UTR_region_length))
 ## end of very fast

#############################################################################
 

######################################################################
 num.rows        <- dim(matr)[1]       

 
 pyrimid         <- c(1,2)
 purin           <- c(3,4)
 XXX             <- new.env()
 XXX$transitions <- vector(,length(biallelic.sites))
 XXX$count       <- 1
 ret.vek         <- integer(num.rows)

  Bial.Mat <- 
       apply(SUBMAT,2,function(x){
       nucs     <- unique(x)
       nuc1     <- x[1]
       nuc2     <- x[2]
       trans1   <- setequal(nucs,pyrimid) # pyrimid
       trans2   <- setequal(nucs,purin)   # purin 
       if(trans1 | trans2){XXX$transitions[XXX$count]<-TRUE}
       XXX$count <- XXX$count + 1 
       nuc1.id   <- nuc1==x
       nuc2.id   <- nuc2==x
       num.nuc1  <- sum(nuc1.id)
       num.nuc2  <- sum(nuc2.id)
       if(num.nuc1<=num.nuc2){
          ret.vek[nuc1.id] <- 1
       }else{
          ret.vek[nuc2.id] <- 1
       }
  return(ret.vek)
  })

# Transitions
  transitions        <- XXX$transitions
  n.transitions      <- sum(transitions)
  n.transversions    <- length(biallelic.sites) - n.transitions
  tt.ratio           <- n.transitions/n.transversions 

return(list(biallelic.matrix=Bial.Mat,biallelic.sites=biallelic.sites,
transitions=transitions,n.valid.sites=NaN,synonymous=rep(NaN,length(biallelic.sites)),
trans.transv.ratio=tt.ratio,Coding_region_length=NaN,Gene_region_length=NaN,Exon_region_length=NaN,
Intron_region_length=NaN,UTR_region_length=NaN))
}

#############################################################################
# END OF FAST

# GFF
# if(is.list(gff)){
# features <- gff
#}



#### INIT
n_site  <- dim(matr)[2]
nsamtot <- dim(matr)[1]

matrix_pol <- integer(nsamtot)



      ### Apply var
      #matrix_pos      <- FALSE
      #mhitbp          <- NaN
      #mis             <- 0
      #matrix_freq     <- NaN
      #matrix_sv       <- NA
      #algsites        <- 0
      #mvbp            <- 0
      #s_sites         <- 0
      #bial_sites      <- 0
      #mis             <- 0
      #gap             <- 0 
      #############

      pyrimid         <- c(1,2)
      purin           <- c(3,4)
      err             <- c(5,6)

GLOBAL <- new.env()
GLOBAL$GLOBALmatrix_pol <- NULL         ### BIALLELIC MATRIX
RETURN                  <-  integer(9)

#count <<- 0

# names(RETURN)    <-  c("matrix_sv","mhitbp","mis","mvbp","algsites","bial_sites")

# 1 : matrix_sv
# 2 : mhitbp
# 3 : mis
# 4 : mvbp
# 5 : algsites
# 6 : bial_sites



### APPLY: ITERATION OVER ALL COLUMNS
erg <- apply(matr,2,function(check){

        
      #cat(count,"\n")
      #count <<- count + 1
 
      fuenfsechs   <- is.element(err,check) # err = c(5,6)
      fuenf        <- fuenfsechs[1]
      sechs        <- fuenfsechs[2]
     
     # keine 5 oder 6
     if(!fuenf & !sechs){
      
       test <- unique(check)
       size <- length(test)

       if(size==1){
          # mono     <- 1
       return(RETURN)
       }

       if(size>2){
         RETURN[5] <- 1 # algsites
         RETURN[2] <- 1 # mhitbp
        return(RETURN)
       }
       
      if(size==2){
       
         RETURN[6] <- 1 # bial_sites
         
         nuc1       <- test[1]
         nuc2       <- test[2]
         nuc12      <- c(nuc1,nuc2)
         # transition/transversion
         trans1        <- setequal(nuc12,pyrimid) # pyrimid
         trans2        <- setequal(nuc12,purin)   # purin        
         ifelse(trans1 | trans2,RETURN[1] <- 1,RETURN[1] <- 0) # matrix_sv
         
         countnuc1  <- sum(check == nuc1)
         countnuc2  <- nsamtot - countnuc1  
     
          if(countnuc1 <= countnuc2){
             RETURN[7]                 <- nuc1
             matrix_pol[check==nuc1]   <- 1 
             GLOBAL$GLOBALmatrix_pol   <- cbind(GLOBAL$GLOBALmatrix_pol,matrix_pol)
            
          }else{
             RETURN[7]                 <- nuc2
             matrix_pol[check==nuc2]   <- 1 
             GLOBAL$GLOBALmatrix_pol   <- cbind(GLOBAL$GLOBALmatrix_pol,matrix_pol)
          }
          
       return(RETURN)    
      } # end size==2
 
    }# end of !fuenfsechs
  

     if(sechs){   
       RETURN[8]   <- 1 # gap
       RETURN[5]   <- 1 # algsites
         
       return(RETURN)
     } # return

    if(fuenf && include.unknown){
    
       mis        <- 1  
       gapids     <- check==5
       check2     <- check[!gapids]
       hh         <- unique(check2)
       size       <- length(hh)
       RETURN[9]  <- 1       

       if(size>2){

         RETURN[5]   <- 1 # algsites
         RETURN[2]   <- 1 # third nucleotide mismatch # mhitbp  
      
       return(RETURN)

       }
       
       if(size==2){ # biallelic
         
         RETURN[4]             <- 1 # mvbp
         RETURN[6]             <- 1 # bial_sites
               
         nuc1       <- hh[1]
         nuc2       <- hh[2]

         nuc12         <- c(nuc1,nuc2)
         # hier transition/transversion
         trans1        <- setequal(nuc12,pyrimid) # pyrimid
         trans2        <- setequal(nuc12,purin)   # purin        
         ifelse(trans1 | trans2,RETURN[1] <- 1,RETURN[1] <- 0) # matrix_sv
         # -----------------------------------------------------------

         countnuc1  <- sum(check2 == nuc1)
         countnuc2  <- sum(check2 == nuc2)
              
          if(countnuc1 <= countnuc2){
             RETURN[7]                 <- nuc1
             matrix_pol[check==nuc1]   <- 1 
          }else{
             RETURN[7]                 <- nuc2
             matrix_pol[check==nuc2]   <- 1 
          }
                     
          matrix_pol[gapids] <- NaN # GAPS
          GLOBAL$GLOBALmatrix_pol   <- cbind(GLOBAL$GLOBALmatrix_pol,matrix_pol)
 
         # biallelic
       } # end size == 2

     return(RETURN)
      
    }# end 
    
    if(fuenf && !include.unknown){
      RETURN[5]          <- 1 # algsites
      RETURN[9]          <- 1
      return(RETURN)
    }
    
 return(RETURN)

}) # End of APPLY

#cat("APPLY is durch \n")

## Function
#####################################
#erg <- as.matrix(erg)
#####################################

# row.names(erg) <- c("matrix_sv","mhitbp","mis","mvbp","algsites","bial_sites")

## BIALLELIC INDICES ######################
matrix_pos <- which(erg[6, ]==1)
bial_sites <- length(matrix_pos)
###########################################

#cat("bial_sites")

### ----------------- Exception ----------------- ####
### if there are no biallelic positions return NA ####
### --------------------------------------------- ####

if(bial_sites==0){
# print("No biallelic positions !")
return(NA)
}
#########################################################

### gaps ################################################
gaps      <- which(erg[8,]==1)
#########################################################

#cat("gaps")

## unknown #############################################
unknown  <- which(erg[9,]==1)
########################################################

#cat("unknown")

## SUM_SAM ##############################################
algsites  <- which(erg[5,]==1)
if(length(algsites)==0){sum_sam  <- rowSums(matr!=5) }else{sum_sam <- rowSums(matr[,-algsites,drop=FALSE]!=5)}
#########################################################

#cat("sum_sam")

### MATRIX_POL (BIALLELIC MATRIX)
matrix_pol           <- GLOBAL$GLOBALmatrix_pol
rm(GLOBAL) # remove environment !
colnames(matrix_pol) <- matrix_pos
#--------------------------------
#################################

#cat("matrix_pol")

###########################################
### Generate codonpositions from matrix_pos
###----------------------------------------

#if(!is.list(gff)){

if(!SNP.DATA){ 

  y <- 1
  matrix_codonpos <- vector(,3*bial_sites)
  
 for(xx in 1:bial_sites){
  x  <- matrix_pos[xx]
  if(x%%3==0){matrix_codonpos[y:(y+2)] <- c(x-2,x-1,x);y <- y+3;next;}
  if(x%%3==1){matrix_codonpos[y:(y+2)] <- c(x,x+1,x+2);y <- y+3;next;}
  if(x%%3==2){matrix_codonpos[y:(y+2)] <- c(x-1,x,x+1);y <- y+3;next;}
 } 

 # matrix_codonpos <- unique(matrix_codonpos) important change !
  
 
 if(matrix_codonpos[length(matrix_codonpos)]>dim(matr)[2]){
    ende            <- length(matrix_codonpos)-3
    matrix_codonpos <- matrix_codonpos[1:ende]
 } # Schmeisse das letzte codon raus ! Kommt nur vor wenn die Original Matrix # #  

}else{matrix_codonpos <- NaN}


#cat("matrix_codon_pos")
##############################################

## GFF
if(is.list(gff)){
   
features <- gff
   
   if(length(features$Gene)>0){
    Gene_region        <- as.vector(unlist(apply(features$Gene,1,function(x){return(x[1]:x[2])})))
    GeneSNPS           <- is.element(matrix_pos,Gene_region)
    Gene_region_length <- length(Gene_region)
    #IntronSNPS    <- matrix_pos[IntronSNPS]
   }else{GeneSNPS<-NaN;Gene_region_length <- 0} 

   if(length(features$Intron)>0){
    Intron_region <- as.vector(unlist(apply(features$Intron,1,function(x){return(x[1]:x[2])})))
    IntronSNPS    <- is.element(matrix_pos,Intron_region)
    Intron_region_length <- length(Intron_region)
    #IntronSNPS    <- matrix_pos[IntronSNPS]
   }else{IntronSNPS<-NaN;Intron_region_length <- 0}

   if(length(features$UTR)>0){	
    UTR_region        <- as.vector(unlist(apply(features$UTR,1,function(x){return(x[1]:x[2])})))
    UTRSNPS           <- is.element(matrix_pos,UTR_region)
    UTR_region_length <- length(UTR_region) 
    #UTRSNPS       <- matrix_pos[UTRSNPS]
   }else{UTRSNPS<-NaN;UTR_region_length <- 0}

   if(length(features$Exon)>0){
    Exon_region        <- as.vector(unlist(apply(features$Exon,1,function(x){return(x[1]:x[2])})))
    ExonSNPS           <- is.element(matrix_pos,Exon_region)
    Exon_region_length <- length(Exon_region)
   }else{ExonSNPS<-NaN;Exon_region_length<-0}
   
   if(length(features$Coding)>0){ 
    Coding_region  <- as.vector(unlist(apply(features$Coding,1,function(x){return(x[1]:x[2])})))
    start.pos      <- as.vector(unlist(apply(features$Coding,1,function(x){return(rep(x[1],x[2]-x[1]+1))})))
    end.pos        <- as.vector(unlist(apply(features$Coding,1,function(x){return(rep(x[2],x[2]-x[1]+1))})))
    # pump reverse strand #FIXME Performance
    JJJ            <- new.env()
    JJJ$count      <- 1
    open(features$rev.strand)
    rev.strand     <-  as.vector(unlist(apply(features$Coding,1,function(x){back <- rep(features$rev.strand[JJJ$count],x[2]-x[1]+1);JJJ$count <- JJJ$count+1;return(back)}))) 
    close(features$rev.strand) 
    #---------------

    ids            <- match(matrix_pos, Coding_region)
    ids            <- ids[!is.na(ids)]
    
    # reverse strand
    rev.strand     <- rev.strand[ids]

 
    start.pos      <- start.pos[ids] 
    end.pos        <- end.pos[ids]
  

    CodingSNPS     <- is.element(matrix_pos,Coding_region)
    Coding_region_length <- length(Coding_region)
    CodingSNPS2     <- matrix_pos[CodingSNPS]
    size            <- length(CodingSNPS2)
   }else{CodingSNPS<-NaN;size<-0;Coding_region_length <- 0}
   
if(size>0 & !SNP.DATA){  # wenn SNPS in den codierenden regionen existieren
   y <- 1 
   matrix_codonpos <- vector(,3*size)
   reverse.codons  <- vector(,3*size)
   rvt             <- c(T,T,T)

 for(xx in 1:size){
  x  <- CodingSNPS2[xx]
  if(rev.strand[xx]){# reverse strand

   if((end.pos[xx]-x)%%3==0){reverse.codons[y:(y+2)] <- rvt;matrix_codonpos[y:(y+2)] <- c(x,x-1,x-2);y <- y+3;next;}
   if((end.pos[xx]-x)%%3==1){reverse.codons[y:(y+2)] <- rvt;matrix_codonpos[y:(y+2)] <- c(x+1,x,x-1);y <- y+3;next;}
   if((end.pos[xx]-x)%%3==2){reverse.codons[y:(y+2)] <- rvt;matrix_codonpos[y:(y+2)] <- c(x+2,x+1,x);y <- y+3;next;}

  }else{# non-reverse strand

   if((x-start.pos[xx])%%3==0){matrix_codonpos[y:(y+2)] <- c(x,x+1,x+2);y <- y+3;next;}
   if((x-start.pos[xx])%%3==1){matrix_codonpos[y:(y+2)] <- c(x-1,x,x+1);y <- y+3;next;}
   if((x-start.pos[xx])%%3==2){matrix_codonpos[y:(y+2)] <- c(x-2,x-1,x);y <- y+3;next;}
 
  }
 
 } 
 
 
  # matrix_codonpos <-  unique(matrix_codonpos) IMPORTANT ! CHECK !
   
 if(matrix_codonpos[length(matrix_codonpos)]>dim(matr)[2]){
    ende            <- length(matrix_codonpos)-3
    matrix_codonpos <- matrix_codonpos[1:ende]
    
 } # Schmeisse das letzte codon raus ! Kommt nur vor wenn die Original Matrix # #  
 
}else{matrix_codonpos <- NaN;CodingSNPS<- NaN;Coding_region_length <- NaN}# if size >0 
}else{CodingSNPS<- NaN;IntronSNPS <- NaN; UTRSNPS <- NaN;ExonSNPS <- NaN;GeneSNPS<-NaN;
Coding_region_length<-NaN;Intron_region_length<-NaN;UTR_region_length<-NaN;Exon_region_length<-NaN;Gene_region_length<-NaN}

# if (bial_sites!=0 & is.list(gff))
# END GFF


### TCGA #################  COUNT THE NUCLEOTIDES PER ROW

nucmat <- matr[,matrix_pos,drop=FALSE]

### Get Nucleotide Substitutions and Mutations #

subst            <- apply(nucmat,2,function(x){return(unique(x[x!=5]))})
colnames(subst)  <- matrix_pos
mutations        <- erg[7,matrix_pos]

#cat("subst")
### -----------------------------#

tcga1  <- rowSums(nucmat==1)
tcga2  <- rowSums(nucmat==2)
tcga3  <- rowSums(nucmat==3)
tcga4  <- rowSums(nucmat==4)
tcga5  <- rowSums(nucmat==5)
tcga6  <- rowSums(nucmat==6)

TCGA           <- cbind(tcga1,tcga2,tcga3,tcga4,tcga5,tcga6)
colnames(TCGA) <- c("T","C","G","A","unknown","GAP")
###############################

#cat("TCGA")

# Syn nonSyn Sites #################

synnonsyn            <- rep(NaN,bial_sites)

if(!SNP.DATA){

 if(!is.na(matrix_codonpos[1]) ){

 # GFF 
 testmatrix           <- matr[,matrix_codonpos,drop=FALSE]

#print(testmatrix)

### change rows for reverse strands
 if(any(reverse.codons)){
 komplement         <- c(4,3,2,1,5,6)   
 testmatrix.reverse <- testmatrix[,reverse.codons,drop=FALSE]
 testmatrix.reverse <- apply(testmatrix.reverse,2,function(x){return(komplement[x])})
 testmatrix[,reverse.codons] <- testmatrix.reverse
 rm(testmatrix.reverse)
 }

 rm(matr)          # -------------------------------------------------------> remove original matrix
 

 colnames(testmatrix) <- matrix_codonpos
 synnonsynL           <- getsyn(testmatrix)
 Codons               <- synnonsynL$Codons
 syn                  <- synnonsynL$synid
 nonsyn               <- synnonsynL$nonsynid
 synid                <- match(syn,matrix_pos)
 synid                <- synid[!is.na(synid)]
 nonsynid             <- match(nonsyn,matrix_pos)
 nonsynid             <- nonsynid[!is.na(nonsynid)]
 synnonsyn[synid]     <- 0 
 synnonsyn[nonsynid]  <- 1 

 }else{syn  <- NaN;nonsyn <- NaN; Codons <- as.list(NaN)}
}else{syn  <- NaN;nonsyn <- NaN; Codons <- as.list(NaN)}
#####################################
 
#cat("synonsyn")

# ALGSITES (SITES WITHOUT VALUES >5 OR MHIT) or ==5 and include.unknown==0)
############################################################################

algsites <- n_site - length(algsites)

############################################################################


##### MHIT (IF THERE IS A >2 NUCLEOTIDE)
mhitbp <- which(erg[2,]==1)
mhit   <- length(mhitbp)
########################
#cat("mhitp")

#### MIS (SUM OF ALL MHITS) # if 5 and includeunknown==1
# totalmis <- sum(erg["mis",])
##############


## MVBP: SITE POSITION -GAP INCLUDED IN THE BIALLELIC MATRIX
mvbp <- which(erg[4,]==1)
###########################

#cat("mvbp")

#### BIALLELIC SUBMATRIX FOR ANALYSIS
erg2 <- erg[ ,matrix_pos,drop=FALSE]
#####################################

## Delete erg ###########################################
rm(erg)
erg <- as.matrix(NaN)
#########################################################


#-------------------------------------#
# VALUES FOR THE BIALLELIC MATRIX !!! #
#-------------------------------------#

# Transition/Transversions ############
# 0: transversions
# 1: transitions
#######################################
matrix_sv      <- erg2[1,]
transitions    <- sum(matrix_sv)
transversions  <- bial_sites - transitions
matrix_sv      <- as.logical(matrix_sv)
#######################################

#cat("transitions")

rm(erg2) ########################### delete erg 2

## matrix_freq #############################
matrix_freq <- colSums(matrix_pol,na.rm=TRUE)
############################################

## UNIC ###################################
unicids      <- which(matrix_freq==1)
unicmatrix   <- matrix_pol[,unicids,drop=FALSE]
unic         <- rowSums(unicmatrix,na.rm=TRUE)
############################################

#cat("singletons")
  #
  
 #print(IntronSNPS)
  
RETURNLISTE <- list(n_site=n_site,transitions=matrix_sv, biallelic.matrix=matrix_pol,
 biallelic.sites=matrix_pos,matrix_codonpos=matrix_codonpos,n.singletons=unic,totalmis=NaN,s_sites=NaN,mvbp=mvbp
,trans.transv.ratio=(transitions/transversions),n.valid.sites=algsites,n.biallelic.sites=bial_sites,
polyallelic.sites=mhitbp,n.nucleotides=sum_sam,biallelic.compositions=TCGA,ROUGH=erg,matrix_freq=matrix_freq,syn=syn,
nonsyn=nonsyn,synonymous=!as.logical(synnonsyn),biallelic.substitutions=subst,minor.alleles=mutations,codons=Codons,
CodingSNPS=CodingSNPS,UTRSNPS=UTRSNPS,IntronSNPS=IntronSNPS,ExonSNPS=ExonSNPS,GeneSNPS=GeneSNPS,Coding_region_length=Coding_region_length,
UTR_region_length=UTR_region_length,Intron_region_length=Intron_region_length,Exon_region_length=Exon_region_length, Gene_region_length=Gene_region_length, sites.with.gaps=gaps,sites.with.unknowns=unknown)


return(RETURNLISTE)

}



