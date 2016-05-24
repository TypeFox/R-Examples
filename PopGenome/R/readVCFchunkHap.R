readVCFchunkHap <- function(vcffile, samplenames=FALSE){


# if(file.exists("SNPRObjects")){unlink("SNPRObjects",recursive=TRUE)}

# Get the sample names !
vcf_handle <- file(vcffile)
open(vcf_handle)

while(1){
 line  <- readLines(vcf_handle,n=1)
 patt  <- substr(line,1,2)
 if(patt!="##"){break}
}

line         <- strsplit(line,"\t")[[1]]
individuals  <- line[10:length(line)]
close(vcf_handle)
#
# End of get individuals

# specify colClasses
spec <- c("NULL","integer","NULL","character",
 "character", rep("NULL",4),rep("character",length(individuals))) 

vcf  <- read.table(vcffile, sep="\t", colClasses=spec)

# Save infos from vcf chunk
pos <- vcf[,1]
ref <- vcf[,2]
alt <- vcf[,3]
ind <- vcf[,4:dim(vcf)[2]]


# FIRST HAPLOTYPE MATRIX !!! #############################################################################

ind_MATRIX <- matrix(,dim(ind)[2],dim(ind)[1])

# fill individual matrix
op <- options(warn = (-1)) # keine Warnmeldungen !
for( xx in 1:dim(ind)[2] ){
  ind_MATRIX[xx,] <- as.integer(substr(ind[,xx], 1, 1)) # 1/1
}
options(op) # Warnmeldungen wieder an !

ind_MATRIX <- ind_MATRIX + 1 


# create ALT MATRIX
ALTlist     <- strsplit(alt,",")
n.val       <- sapply(ALTlist,function(x){return(length(x))})
 

alt_MATRIX <- matrix(,length(pos),max(n.val))
for(xx in 1:length(ALTlist)){
alt_MATRIX[xx,1:length(ALTlist[[xx]])] <- ALTlist[[xx]] 
}


# MAP Matrix
MAP_MATRIX <- cbind(ref,alt_MATRIX)


Code_matrix <- matrix(,dim(ind_MATRIX)[1], length(pos)) 
for(xx in 1:dim(MAP_MATRIX)[1]){
 Code_matrix[,xx] <- MAP_MATRIX[xx,][ind_MATRIX[,xx]]
}


Code_matrix <- .Call("code_nucs",Code_matrix)   
# ----------------------------------------------

# Set the rownames
if(samplenames[1]!=FALSE){
   rownames(Code_matrix) <- samplenames
}else{
   rownames(Code_matrix) <- individuals
}

o_b_j_sub    <- list(matrix=Code_matrix,reference=NaN,positions=pos)

return(o_b_j_sub)
}

