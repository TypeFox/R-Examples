PopGenread <- function(filepath,format) {

if(format=="VCF"){
 res  <- myReadVCF(filepath)
 mat  <- res$matrix
 ref  <- res$reference
 pos  <- res$positions
 return(list(matrix=mat,reference=ref,positions=pos)) 
}

############################################################
#if(format=="VCF"){
# res  <- readVCFchunk(filepath)
# mat  <- res$matrix
# ref  <- res$reference
# pos  <- res$positions
# return(list(matrix=mat,reference=ref,positions=pos)) 
#}

#if(format=="VCFhap"){
# res  <- readVCFchunkHap(filepath)
# mat  <- res$matrix
# ref  <- res$reference
# pos  <- res$positions
# return(list(matrix=mat,reference=ref,positions=pos)) 
#}

#if(format=="VCFtri"){
# res  <- readVCFchunk_tri(filepath)
# mat  <- res$matrix
# ref  <- res$reference
# pos  <- res$positions
# return(list(matrix=mat,reference=ref,positions=pos)) 
#}

#if(format=="VCFtet"){
# res  <- readVCFchunk_tet(filepath)
# mat  <- res$matrix
# ref  <- res$reference
# pos  <- res$positions
# return(list(matrix=mat,reference=ref,positions=pos)) 
#}
###############################################################

if(format=="HapMap"){
 res  <- parse_HapMap(filepath)
 mat  <- res$matrix
 ref  <- res$reference
 pos  <- res$positions
 return(list(matrix=mat,reference=ref,positions=pos)) 
}

if(format=="RData"){
 XXX   <- load(filepath)
 mmm   <- get(XXX[1])
 return(list(matrix=mmm$matrix,reference=NaN,positions=mmm$positions))
}


if(format=="nexus"){

 matrix     <- my_read.nexus(filepath)
 nn         <- rownames(matrix)

 number     <- c(1,1,1,1,2,2,3,3,4,4,5,5,5,6)
 nuc        <- c("T","t","U","u","C","c","G","g","A","a","N","n","?","-")


  matrix <- apply(matrix,1,function(x){return(as.integer(number[match(x,nuc)]))})
  matrix <- t(matrix)

  matrix[is.na(matrix)] <- 5
  rownames(matrix)      <- nn
  attr(matrix,"path")   <- filepath
return(list(matrix=matrix,reference=NaN,positions=NaN))
}


   gen            <-  .Call("readdna",filepath)
   #gen            <-  .Call("my_read_fasta",filepath)
   rownames(gen)  <-  gsub(" ","",rownames(gen))
   rownames(gen)  <-  gsub("\r","",rownames(gen))
   if(is.null(gen)){gen <- NaN}  

#----ape
#  gen     <- read.dna(filepath,format,as.character=TRUE) # ape package
#  r.names <- rownames(gen)
#  gen     <- .Call("code_nucs",gen)
#  rownames(gen) <- r.names
#----ape  

  return(list(matrix=gen,reference=NaN,positions=NaN))

#############################################################
  
 
}
