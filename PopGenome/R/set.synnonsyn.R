setGeneric("set.synnonsyn", function(object,ref.chr=FALSE,save.codons=FALSE) standardGeneric("set.synnonsyn"))
 setMethod("set.synnonsyn", "GENOME",
 function(object,ref.chr,save.codons){


if(ref.chr[1]==FALSE){
stop("Please verify the reference sequence")
}

if(object@gff.info == FALSE) {
stop("No GFF file was read in !")
}


for (xyz in 1:length(ref.chr)){
# erstmal nur fuer ein Chunk

Coding.matrix            <- object@region.data@Coding.matrix2[[xyz]][,] # weil ff object, 2 because (fitting GFF)
biallelic.sites2         <- object@region.data@biallelic.sites2[[xyz]]  # with respect to the refernce positions
biallelic.sites          <- object@region.data@biallelic.sites[[xyz]]
START                    <- object@region.data@reading.frame[[xyz]][,] #  weil ff object
REV                      <- object@region.data@rev.strand[[xyz]][,]    #  reverse strand information #---
REV                      <- as.logical(REV) 
CodingSNPS               <- object@region.data@CodingSNPS[[xyz]]

if(sum(CodingSNPS)==0){
warning("No coding SNPs in this region !")
print(object@region.names[xyz])
next
}

#print(REV)
if(any(!REV)){
pos.strand.shift         <- START[!REV,3]

#zero			 <- pos.strand.shift==0
#one 			 <- pos.strand.shift==1
#two 			 <- pos.strand.shift==2	
#pos.strand.shift[zero]   <- 0
#pos.strand.shift[one]    <- 2
#pos.strand.shift[two]    <- 1

}else{
pos.strand.shift         <- 0
}

if(any(REV)){
rev.strand.shift         <- START[REV,3] 

#zero			 <- rev.strand.shift==0 
#one 			 <- rev.strand.shift==1 
#two 			 <- rev.strand.shift==2
#rev.strand.shift[zero]   <-  0
#rev.strand.shift[one]    <-  2
#rev.strand.shift[two]    <-  1
}else{
rev.strand.shift         <- 0
}

START                    <- START[,1]
START[!REV]              <- START[!REV] + pos.strand.shift

START[REV]               <- object@region.data@reading.frame[[xyz]][REV,2] - rev.strand.shift  

#START[REV]  - rev.strand.shift


#START                    <- START[,1] + START[,2]

Coding.matrix            <- Coding.matrix


#print(START)


# ---
# unique
# wenn regionen doppelt
#ids                     <- !duplicated(Coding.matrix[,1])
#START                   <- START[ids]
#Coding.matrix           <- Coding.matrix[ids,,drop=FALSE]

# define an evironment
synGLOBAL <- new.env() 

# Create Region and save size of region 
synGLOBAL$SIZE  <- numeric(dim(Coding.matrix)[1])
synGLOBAL$count <- 1

erg  <- apply(Coding.matrix,1,function(xx){

 region                            <- xx[1]:xx[2] 
 synGLOBAL$SIZE[synGLOBAL$count]   <- length(region)
 synGLOBAL$count                   <- synGLOBAL$count + 1 
 return(region)

})

# what are the real positions ?
 erg         <- unlist(erg)

 bial.pos    <- match(erg,biallelic.sites2) #.Call("my_match_C",erg,biallelic.sites2)
 bial.pos[bial.pos==-1] <- NaN
 #return(bial.pos)
 bial.pos    <- biallelic.sites[bial.pos]

#print(bial.pos)
#print(biallelic.sites2)

# Create Start Vector
 synGLOBAL$count <- 1
 vec <- sapply(START,function(x){
      gg              <- rep(x,synGLOBAL$SIZE[synGLOBAL$count])       
      synGLOBAL$count <- synGLOBAL$count + 1

 return(gg)
 })

# Create REV Vector #----
 synGLOBAL$count <- 1
 vec_rev <- sapply(REV,function(x){
      gg              <- rep(x,synGLOBAL$SIZE[synGLOBAL$count])       
      synGLOBAL$count <- synGLOBAL$count + 1

 return(gg)
 })


 START.vec   <- unlist(vec)
 REV.vec     <- unlist(vec_rev) #---

#print(length(REV.vec))
#print(length(START.vec))
#print(length(bial.pos))

 cod.pos              <- (bial.pos - START.vec)%%3
 # in case of reverse strand 
 cod.pos[REV.vec]     <- (START.vec[REV.vec]-bial.pos[REV.vec])%%3
 #FIXME


  

 # Delete NaNs 
 cod.pos     <- cod.pos[!is.na(bial.pos)]
 rev.pos     <- REV.vec[!is.na(bial.pos)] #----
 bial.pos    <- bial.pos[!is.na(bial.pos)]

 
 #print(START.vec[300])

# print(length(bial.pos))
# print(length(cod.pos))
# print(length(rev.pos))

 ids         <- !duplicated(bial.pos)
 cod.pos     <- cod.pos[ids]
 bial.pos    <- bial.pos[ids]
 rev.pos     <- rev.pos[ids] #----

# print(bial.pos)
# print(length(cod.pos))
# print(length(rev.pos))

 # print(rev.pos[1215])
 # print(bial.pos[1215])
 # print(cod.pos[1215])
 #print(rev.pos)

# Create the codons
# bial.pos and cod.pos

codons <- matrix(,length(cod.pos),3)

for (xx in 1:length(cod.pos)){
   
    if(rev.pos[xx]){# reverse strand #FIXME
     if(cod.pos[xx]==0){codons[xx,]=c(bial.pos[xx],bial.pos[xx]-1,bial.pos[xx]-2);next}
     if(cod.pos[xx]==1){codons[xx,]=c(bial.pos[xx]+1,bial.pos[xx],bial.pos[xx]-1);next}
     if(cod.pos[xx]==2){codons[xx,]=c(bial.pos[xx]+2,bial.pos[xx]+1,bial.pos[xx]);next}
    }else{
     if(cod.pos[xx]==0){codons[xx,]=c(bial.pos[xx],bial.pos[xx]+1,bial.pos[xx]+2);next}
     if(cod.pos[xx]==1){codons[xx,]=c(bial.pos[xx]-1,bial.pos[xx],bial.pos[xx]+1);next}
     if(cod.pos[xx]==2){codons[xx,]=c(bial.pos[xx]-2,bial.pos[xx]-1,bial.pos[xx]);next}
     }

}
#print(cod.pos)
#print(codons)

## Reading the reference chromosome
file.info <- .Call("get_dim_fasta",ref.chr[xyz])

gc()
#print(file.info)

CHR       <- .Call("get_ind_fasta",ref.chr,1,file.info[[1]][2])

#print(CHR[1:10])

# Create codons with nucleotides
Nuc.codons    <- CHR[codons]
Nuc.codons    <- matrix(Nuc.codons,ncol=3)

ALT           <- Nuc.codons
REF           <- Nuc.codons
Subst         <- object@region.data@biallelic.substitutions[[xyz]]
minor         <- Subst[1,CodingSNPS]
major         <- Subst[2,CodingSNPS]

komplement <- c(4,3,2,1,5)

#print(Nuc.codons)

for(xx in 1: dim(Nuc.codons)[1]){
 

 if(rev.pos[xx]){

  #convert to komplement nucleotides
  REF[xx,]  <- komplement[REF[xx,]]
  ALT[xx,]  <- komplement[ALT[xx,]]
  minor[xx] <- komplement[minor[xx]]
  major[xx] <- komplement[major[xx]]
  ###########
  
  if(cod.pos[xx]==0){REF[xx,1] <- minor[xx];ALT[xx,1]<-major[xx];next}
  if(cod.pos[xx]==1){REF[xx,2] <- minor[xx];ALT[xx,2]<-major[xx];next}
  if(cod.pos[xx]==2){REF[xx,3] <- minor[xx];ALT[xx,3]<-major[xx];next}
 }else{
  if(cod.pos[xx]==0){REF[xx,1] <- minor[xx];ALT[xx,1]<-major[xx];next}
  if(cod.pos[xx]==1){REF[xx,2] <- minor[xx];ALT[xx,2]<-major[xx];next}
  if(cod.pos[xx]==2){REF[xx,3] <- minor[xx];ALT[xx,3]<-major[xx];next}
 }

} 

#print(REF)
#print(ALT)

# Coding Codons ...



ALT <- codonise64(ALT)
REF <- codonise64(REF)

if(save.codons){
saveALTREF <- cbind(REF,ALT)
}

CC  <- codontable()

ALT <- CC$Protein[1,ALT]
REF <- CC$Protein[1,REF]

CHECK <- cbind(ALT,REF)



erg <- apply(CHECK,1,function(x){return(length(unique(x)))})
erg[erg==2] <- 0 #nonsyn
erg[erg==1] <- 1 #syn

# Change object of class GENOME
change <- object@region.data
change@synonymous[[xyz]][CodingSNPS] <- erg

### save codons
if(save.codons){
n.coding.snps     <- sum(CodingSNPS)
codonlist  <- vector("list",n.coding.snps)
count <- 1

for(vv in 1:n.coding.snps){
    codonlist[[vv]] <- saveALTREF[vv,]
}
change@codons[[xyz]] <- codonlist
}
#######################

object@region.data <- change 
}# End Iteration over chunks or chromosomes

return(object)

})
