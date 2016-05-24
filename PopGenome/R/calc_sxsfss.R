calc_sxsfss <- function(matrix_pol,populations,outgroup=FALSE,data){

if(outgroup[1]==FALSE){outgroup <- NULL}
#### if only one population is defined 
if(length(populations)== 1 & (length(outgroup)==0|length(outgroup)==0)){return(list(NaN))}
#----------------------------------------------------------------------------------
if(length(outgroup)==0){outgroup<-NULL}

npops <-  length(populations)
erg   <-  apply(matrix_pol,2,ispolmisX) # look at the whole matrix 
rownames(erg) <- c("bial?","anc")

id    <-  which(erg["bial?",]==1) # 1 is polymorphic

### only polymorphic sites #### its not necessary because all columns in matrix_pol are polymorph
matrix_pol <- matrix_pol[,id,drop=FALSE]
colnames(matrix_pol) <- id
#-----------------------------#

# ------------------------------
# Check the outgroup 
# Check if monomorph or not 
# ------------------------------
if(length(outgroup)>0){ # Outgroup is defined

 if(is.matrix(matrix_pol[outgroup,])){ outgr <- apply(matrix_pol[outgroup,],2,ispolmisX)}
 if(is.vector(matrix_pol[outgroup,])){ outgr <- sapply(matrix_pol[outgroup,],ispolmisX)}
 rownames(outgr) <- c("bial?","anc")

}

################################
##  Check the populations
##  Check if monomorphic or not.
##  polqa ######################

popmatrix <- NULL
for(xx in 1:npops){
    pop        <- apply(matrix_pol[populations[[xx]],,drop=FALSE],2,ispolmisX)
    popmatrix  <- rbind(popmatrix,pop)
}# End of for

  # INIT rownames-#
  popnames  <- paste("pop",1:npops)
  ancnames  <- paste("anc",1:npops)
  NAMES     <- as.vector(mapply(union,popnames,ancnames))
  rownames(popmatrix) <- NAMES
  colnames(popmatrix) <- id
  # ------------- #


ancmatrix <- popmatrix[ancnames,,drop=FALSE]
popmatrix <- popmatrix[popnames,,drop=FALSE]


### polqb ##### ------------------------------------------------------------ # polqb #
if(npops>1){ # more than one population  --> calculate poplqb
ohnepopmatrix  <- NULL


for(xx in 1:npops){

  without       <- populations[[xx]]
  pops          <- unlist(populations)
  Xid           <- match(without,pops)
  ohnepop       <- pops[-Xid]

  ohnexx        <- matrix_pol[ohnepop,,drop=FALSE]            # Betrachte die restlichen Pops ausser xx
  poprest       <- apply(ohnexx,2,ispolmisX)
  ohnepopmatrix <- rbind(ohnepopmatrix,poprest)

}
rownames(ohnepopmatrix) <- NAMES
colnames(ohnepopmatrix) <- id

ancmatrixrest <- ohnepopmatrix[ancnames,,drop=FALSE]
popmatrixrest <- ohnepopmatrix[popnames,,drop=FALSE]

}# End if npops > 1 # ------------------------------------------------------ # 

#init ---------------------------
comp1 <- vector("list",npops)
comp2 <- vector("list",npops)
comp3 <- vector("list",npops)
comp4 <- vector("list",npops)
comp5 <- vector("list",npops)
#-------------------------------

 # attr(SX_SXF_SF_SS,"SX :") <-   "I am polymorph rest is monomorph "
 # attr(SX_SXF_SF_SS,"SXF:") <-   "I am monomorph rest is polymorph "
 # attr(SX_SXF_SF_SS,"SF:")  <-   "I am mono rest is mono with same mono value"
 # attr(SX_SXF_SF_SS,"SS:")  <-   "I am mono rest is mono with different mono value" 
 # attr(SX_SXF_SF_SS,"SXX:") <-   "I am poly rest is poly" 

#-------------- Outgroup is defined ---------------------------------------- #
if(length(outgroup)>0){
 polc      <- outgr[1,]  # polymorphic 1 or not 0
 polc_anc  <- outgr[2,]  # outgroup mono value
}
#--------------------------------------------------------------------------- #

if(npops > 1 & length(outgroup)>0){ # mehr als eine Population und Outgroup definiert

for(xx in 1:npops){

   polqa      <- popmatrix[xx,]
   polqa_anc  <- ancmatrix[xx,]
   polqb      <- popmatrixrest[xx,]
   polqb_anc  <- ancmatrixrest[xx,]

   # Wenn keine outgroup definiert polc = NaN
   comp1[[xx]]  <- as.numeric(names(which(polqa==1 & polqb==0 & polc==0 & (polqb_anc==polc_anc))))  # xx is polymorph rest is (monomorph) and identical with the OUTGROUP !!!
   comp2[[xx]]  <- as.numeric(names(which(polqa==1 & polqb==0 & polc==0 &  polqb_anc!=polc_anc)))   # xx is polymorph rest is monomorph, but not identical with the outgroup
   comp3[[xx]]  <- as.numeric(names(which(polqa==0 & polqb==0 & polc==0 &  polqb_anc==polc_anc)))   # xx is monomorph rest is monomorph and identical with the outgroup !!!
   comp4[[xx]]  <- as.numeric(names(which(polqa==1 & polqb==1 & polc==0))) # xx is polymorph , the rest is polymorph, outgroup is monomorph

   comp1[[xx]] <- data$biallelic.sites[comp1[[xx]]]
   comp2[[xx]] <- data$biallelic.sites[comp2[[xx]]]
   comp3[[xx]] <- data$biallelic.sites[comp3[[xx]]]
   comp4[[xx]] <- data$biallelic.sites[comp4[[xx]]]
}
} # End if if(npops > 1 & length(outgroup)>0)

if(npops==1 & length(outgroup)>0){ # mehr als eine Population und Outgroup ist nicht definiert
   
   polqa      <- popmatrix
   polqa_anc  <- ancmatrix
   
   comp1[[xx]]  <- as.numeric(names(which(polqa==1 & polc==0 )))  # xx is polymorph Outgroup is (monomorph) 
   comp2[[xx]]  <- as.numeric(names(which(polqa==0 & polc==1)))   # xx is monomorph Outgroup is polymorph
   comp3[[xx]]  <- as.numeric(names(which(polqa==0 & polc==0 & polqa_anc==polc_anc)))  # xx is monomorph Outgroup is monomorph and have the same monomorph value
   comp4[[xx]]  <- as.numeric(names(which(polqa==0 & polc==0 & polqa_anc!=polc_anc)))  # xx is monomorph Outgroup is monomorph but not with the same monomorph value 

   comp1[[xx]] <- data$biallelic.sites[comp1[[xx]]]
   comp2[[xx]] <- data$biallelic.sites[comp2[[xx]]]
   comp3[[xx]] <- data$biallelic.sites[comp3[[xx]]]
   comp4[[xx]] <- data$biallelic.sites[comp4[[xx]]]
   
} # mehr als eine Population, aber keine Outgroup definiert

if(npops > 1 & length(outgroup)==0){ # mehr als eine Population, aber keine Outgroup ist definiert

  for(xx in 1:npops){
   
   polqa      <- popmatrix[xx,]
   polqa_anc  <- ancmatrix[xx,]
   polqb      <- popmatrixrest[xx,]
   polqb_anc  <- ancmatrixrest[xx,]
   
   comp1[[xx]]  <- as.numeric(names(which(polqa==1 & polqb==0 )))  # xx is polymorph rest is (monomorph) 
   comp2[[xx]]  <- as.numeric(names(which(polqa==0 & polqb==1)))   # xx is monomorph rest is polymorph
   comp3[[xx]]  <- as.numeric(names(which(polqa==0 & polqb==0 & polqa_anc==polqb_anc))) # xx is monomorph rest is monomorph and have the same monomorph value
   comp4[[xx]]  <- as.numeric(names(which(polqa==0 & polqb==0 & polqa_anc!=polqb_anc))) # xx is monomorph rest is monomorph but not with the same monomorph value 
   
   # comp5[[xx]]  <- as.numeric(names(which(polqa==1 & polqb==1))) #  xx is polymorph rest is polymorph

   comp1[[xx]] <- data$biallelic.sites[comp1[[xx]]]
   comp2[[xx]] <- data$biallelic.sites[comp2[[xx]]]
   comp3[[xx]] <- data$biallelic.sites[comp3[[xx]]]
   comp4[[xx]] <- data$biallelic.sites[comp4[[xx]]]
   
   # comp5[[xx]] <- data$biallelic.sites[comp5[[xx]]]
  
  }

} # mehr als eine Population, aber keine Outgroup definiert

if(length(outgroup>0) & npops>1){
   
   compout <- vector("list",3)
   
   # ---- OUTGROUP ----------------------------------
   # Look at the last population
    
    ONE   <- polc==1 &  polqa==0 & polqb==0   &  (polqb_anc==polqa_anc | npops==1)      # Sxo    ---    Outgroup:polymorph, Populationen: alle identisch
    TWO   <- polc==1 &  polqa==0 & polqb==0   &   polqb_anc!=polqa_anc                  # Sanco  ---    Outgroup:polymorph, Pops:monomorph, popxxancester untersch..
    THREE <- polc==1 & (polqa!=0 | polqb!=0)  &  (polqa!=-1 & polqb!=-1)                # Sanco  ---
    FOUR  <- polc==1 & (polqa!=0 | polqb!=0)  &  (polqa==-1 | polqb==-1) & (polqa==0 | polqb==0) # sxo---- Outgroup:polymorph, 
    FIVE  <- polc==1 & (polqa!=0 | polqb!=0)  &  (polqa==-1 | polqb==-1) & (polqa!=0 & polqb!=0) # sanco
    SIX   <- polc==0 &  polqa==0 & polqb==0   &   polqb_anc==polqa_anc   & polqa_anc!=polc_anc # Sfo  --- Outgroup:monomorph, Pops: identisch ! aber unterschiedlich zur Outgroup
      
    compout[[1]] <- as.numeric(names(which( ONE | FOUR ))) #sxo
    compout[[2]] <- as.numeric(names(which( TWO | THREE | FIVE ))) #sanco
    compout[[3]] <- as.numeric(names(which( SIX ))) #Sfo 
    
    compout[[1]] <- data$biallelic.sites[compout[[1]]]
    compout[[2]] <- data$biallelic.sites[compout[[2]]]
    compout[[3]] <- data$biallelic.sites[compout[[3]]]
    
    compout <- as.matrix(compout)    

}# End if Outroup is defined

# --------- RETURN VALUES -------------------------------------- #

if(npops>1 & length(outgroup)>0){
 SX      <- as.matrix(comp1)
 SXF     <- as.matrix(comp2)
 SF      <- as.matrix(comp3)
 SS      <- as.matrix(comp4)
 
 SX_SXF_SF_SS <- cbind(SX,SXF,SF,SS)
 rownames(SX_SXF_SF_SS) <- popnames
 colnames(SX_SXF_SF_SS) <- c("SX","SXF","SF","SS")
 attr(SX_SXF_SF_SS,"POPULATIONS:") <- ""
 OUT           <- compout
 colnames(OUT) <- "Outgroup"
 rownames(OUT) <- c("Sxo","Sanco","Sfo") 
 attr(OUT,"OUTGROUP:") <- "YES"
 return(list(POP=SX_SXF_SF_SS,OUT=OUT))
}

if(npops>1 & length(outgroup)==0){
 SX    <- as.matrix(comp1)
 SXF   <- as.matrix(comp2)
 SF    <- as.matrix(comp3)
 SS    <- as.matrix(comp4)
 # SXX <- as.matrix(comp5)

 SX_SXF_SF_SS <- cbind(SX,SXF,SF,SS)
 rownames(SX_SXF_SF_SS) <- popnames
 colnames(SX_SXF_SF_SS) <- c("SX","SXF","SF","SS")
 attr(SX_SXF_SF_SS,"POP:") <- " >1"
 attr(SX_SXF_SF_SS,"OUTGROUP") <- "NO" 
 attr(SX_SXF_SF_SS,"SX :") <-  "I am polymorph rest is monomorph "
 attr(SX_SXF_SF_SS,"SXF:") <-  "I am monomorph rest is polymorph "
 attr(SX_SXF_SF_SS,"SF:") <-   "I am mono rest is mono with same mono value"
 attr(SX_SXF_SF_SS,"SS:") <-   "I am mono rest is mono with different mono value" 
 return(list(POP=SX_SXF_SF_SS))
}

if(npops==1 & length(outgroup)>= 1){
 SX    <- as.matrix(comp1)
 SXF   <- as.matrix(comp2)
 SF    <- as.matrix(comp3)
 SS    <- as.matrix(comp4)
 SX_SXF_SF_SS <- cbind(SX,SXF,SF,SS)
 rownames(SX_SXF_SF_SS) <- popnames
 colnames(SX_SXF_SF_SS) <- c("SX","SXF","SF","SS")
 attr(SX_SXF_SF_SS,"POP:") <- " == 1"
 attr(SX_SXF_SF_SS,"OUTGROUP") <- "YES" 
 attr(SX_SXF_SF_SS,"SX :") <-  "I am polymorph Outgroup is monomorph "
 attr(SX_SXF_SF_SS,"SXF:") <-  "I am monomorph Outgroup is polymorph "
 attr(SX_SXF_SF_SS,"SF:") <-   "I am mono Outgroup is mono with same mono value"
 attr(SX_SXF_SF_SS,"SS:") <-   "I am mono Outgroup is mono with different mono value" 
 
 return(list(POP=SX_SXF_SF_SS))
}

}# End of Function

######################################################################
### FUNCTION: ISPOLMISX
### Already implemented in a another function # here with other value
######################################################################

ispolmisX <- function(vek){
  gapids <- !is.na(vek)
  vek    <- vek[gapids]
  ss     <- unique(vek)

  if(length(ss)==1) {return(c(0,as.numeric(ss)))}   #  monomorph
  if(length(ss)==2) {return(c(1,NaN))}              #  biallelic
  if(length(ss)==0) {return(c(-1,NaN))}             #  all are gaps
}
