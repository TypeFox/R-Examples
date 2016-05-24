readVCFchunk <- function(vcffile, samplenames=FALSE){

### DIPLOID DATA

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


# Create Coding 
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

HAPMATRIX1 <- Code_matrix
rm(Code_matrix)

### END OF FIRST HAPLOTYPE MATRIX !!! ##############################################################


### SECOND HAPLOTYPE MATRIX ########################################################################
ind_MATRIX <- matrix(,dim(ind)[2],dim(ind)[1])
# fill individual matrix
op <- options(warn = (-1)) # keine Warnmeldungen !
for( xx in 1:dim(ind)[2] ){
  ind_MATRIX[xx,] <- as.integer(substr(ind[,xx], 3, 3)) # 1/1
}
options(op) # Warnmeldungen wieder an !

ind_MATRIX <- ind_MATRIX + 1 

Code_matrix <- matrix(,dim(ind_MATRIX)[1], length(pos)) 
for(xx in 1:dim(MAP_MATRIX)[1]){
 Code_matrix[,xx] <- MAP_MATRIX[xx,][ind_MATRIX[,xx]]
}

Code_matrix <- .Call("code_nucs",Code_matrix)   

gc()
# ----------------------------------------------

# Set the rownames
if(samplenames[1]!=FALSE){
   rownames(Code_matrix) <- paste(samplenames,".2",sep="")
}else{
   rownames(Code_matrix) <- paste(individuals,".2",sep="")
}

HAPMATRIX2 <- Code_matrix

### END OF SECOND HAPLOTYPE MATRIX ################################################################


Code_matrix            <- rbind(HAPMATRIX1,HAPMATRIX2)

rm(HAPMATRIX1)
rm(HAPMATRIX2)
gc()

# modify row Order
NN           <- rownames(Code_matrix)
names(NN)    <- 1:length(NN)
NN           <- sort(NN)
ids          <- as.integer(names(NN))
Code_matrix  <- Code_matrix[ids,,drop=FALSE]
names(NN)    <- NULL
rownames(Code_matrix) <- NN
# ---------------

o_b_j_sub    <- list(matrix=Code_matrix,reference=NaN,positions=pos)

return(o_b_j_sub)
}
# END DIPLOID DATA


# Outsource data # ---------------------------------------
#  dir.create("SNPRObjects")
 
#  o_b_j_sub    <- list(matrix=Code_matrix,reference=NaN,positions=pos)
#  save(o_b_j_sub,file=file.path(getwd(),"SNPRObjects",paste("chr",sep="")))

#---------------------------------------------------------

#  path     <- file.path(getwd(),"SNPRObjects")
#  res      <- readData(path,format="RData",big.data=TRUE,SNP.DATA=TRUE)

#unlink("SNPRObjects", recursive=TRUE)
#return(res)

#}
###############################################################################################################

readVCFchunk_tet <- function(vcffile, samplenames=FALSE){

### Tetraploid DATA

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


# Create Coding 
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

HAPMATRIX1 <- Code_matrix
rm(Code_matrix)

### END OF FIRST HAPLOTYPE MATRIX !!! ##############################################################


### SECOND HAPLOTYPE MATRIX ########################################################################
ind_MATRIX <- matrix(,dim(ind)[2],dim(ind)[1])
# fill individual matrix
op <- options(warn = (-1)) # keine Warnmeldungen !
for( xx in 1:dim(ind)[2] ){
  ind_MATRIX[xx,] <- as.integer(substr(ind[,xx], 3, 3)) # 1/1
}
options(op) # Warnmeldungen wieder an !

ind_MATRIX <- ind_MATRIX + 1 

Code_matrix <- matrix(,dim(ind_MATRIX)[1], length(pos)) 
for(xx in 1:dim(MAP_MATRIX)[1]){
 Code_matrix[,xx] <- MAP_MATRIX[xx,][ind_MATRIX[,xx]]
}

Code_matrix <- .Call("code_nucs",Code_matrix)   

gc()
# ----------------------------------------------

# Set the rownames
if(samplenames[1]!=FALSE){
   rownames(Code_matrix) <- paste(samplenames,".2",sep="")
}else{
   rownames(Code_matrix) <- paste(individuals,".2",sep="")
}

HAPMATRIX2 <- Code_matrix

### END OF SECOND HAPLOTYPE MATRIX ################################################################

### THIRD HAPLOTYPE MATRIX ########################################################################
ind_MATRIX <- matrix(,dim(ind)[2],dim(ind)[1])
# fill individual matrix
op <- options(warn = (-1)) # keine Warnmeldungen !
for( xx in 1:dim(ind)[2] ){
  ind_MATRIX[xx,] <- as.integer(substr(ind[,xx], 5, 5)) # 1/1
}
options(op) # Warnmeldungen wieder an !

ind_MATRIX <- ind_MATRIX + 1 

Code_matrix <- matrix(,dim(ind_MATRIX)[1], length(pos)) 
for(xx in 1:dim(MAP_MATRIX)[1]){
 Code_matrix[,xx] <- MAP_MATRIX[xx,][ind_MATRIX[,xx]]
}

Code_matrix <- .Call("code_nucs",Code_matrix)   

gc()
# ----------------------------------------------

# Set the rownames
if(samplenames[1]!=FALSE){
   rownames(Code_matrix) <- paste(samplenames,".3",sep="")
}else{
   rownames(Code_matrix) <- paste(individuals,".3",sep="")
}

HAPMATRIX3 <- Code_matrix

### END OF THIRD HAPLOTYPE MATRIX ################################################################

### FOURTH HAPLOTYPE MATRIX ########################################################################
ind_MATRIX <- matrix(,dim(ind)[2],dim(ind)[1])
# fill individual matrix
op <- options(warn = (-1)) # keine Warnmeldungen !
for( xx in 1:dim(ind)[2] ){
  ind_MATRIX[xx,] <- as.integer(substr(ind[,xx], 7, 7)) # 1/1
}
options(op) # Warnmeldungen wieder an !

ind_MATRIX <- ind_MATRIX + 1 

Code_matrix <- matrix(,dim(ind_MATRIX)[1], length(pos)) 
for(xx in 1:dim(MAP_MATRIX)[1]){
 Code_matrix[,xx] <- MAP_MATRIX[xx,][ind_MATRIX[,xx]]
}

Code_matrix <- .Call("code_nucs",Code_matrix)   

gc()
# ----------------------------------------------

# Set the rownames
if(samplenames[1]!=FALSE){
   rownames(Code_matrix) <- paste(samplenames,".4",sep="")
}else{
   rownames(Code_matrix) <- paste(individuals,".4",sep="")
}

HAPMATRIX4 <- Code_matrix

### END OF FOURTH HAPLOTYPE MATRIX ################################################################

Code_matrix            <- rbind(HAPMATRIX1,HAPMATRIX2,HAPMATRIX3,HAPMATRIX4)

rm(HAPMATRIX1)
rm(HAPMATRIX2)
rm(HAPMATRIX3)
rm(HAPMATRIX4)
gc()

# modify row Order
NN           <- rownames(Code_matrix)
names(NN)    <- 1:length(NN)
NN           <- sort(NN)
ids          <- as.integer(names(NN))
Code_matrix  <- Code_matrix[ids,,drop=FALSE]
names(NN)    <- NULL
rownames(Code_matrix) <- NN
# ---------------

o_b_j_sub    <- list(matrix=Code_matrix,reference=NaN,positions=pos)

return(o_b_j_sub)
}

# Outsource data # ---------------------------------------
#  dir.create("SNPRObjects")
 
#  o_b_j_sub    <- list(matrix=Code_matrix,reference=NaN,positions=pos)
#  save(o_b_j_sub,file=file.path(getwd(),"SNPRObjects",paste("chr",sep="")))

#---------------------------------------------------------

#  path     <- file.path(getwd(),"SNPRObjects")
#  res      <- readData(path,format="RData",big.data=TRUE,SNP.DATA=TRUE)

#unlink("SNPRObjects", recursive=TRUE)
#return(res)

#}
#################################################################################################END OF TETRAPLOID DATA

readVCFchunk_tri <- function(vcffile, samplenames=FALSE){

### Tetraploid DATA

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


# Create Coding 
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

HAPMATRIX1 <- Code_matrix
rm(Code_matrix)

### END OF FIRST HAPLOTYPE MATRIX !!! ##############################################################


### SECOND HAPLOTYPE MATRIX ########################################################################
ind_MATRIX <- matrix(,dim(ind)[2],dim(ind)[1])
# fill individual matrix
op <- options(warn = (-1)) # keine Warnmeldungen !
for( xx in 1:dim(ind)[2] ){
  ind_MATRIX[xx,] <- as.integer(substr(ind[,xx], 3, 3)) # 1/1
}
options(op) # Warnmeldungen wieder an !

ind_MATRIX <- ind_MATRIX + 1 

Code_matrix <- matrix(,dim(ind_MATRIX)[1], length(pos)) 
for(xx in 1:dim(MAP_MATRIX)[1]){
 Code_matrix[,xx] <- MAP_MATRIX[xx,][ind_MATRIX[,xx]]
}

Code_matrix <- .Call("code_nucs",Code_matrix)   

gc()
# ----------------------------------------------

# Set the rownames
if(samplenames[1]!=FALSE){
   rownames(Code_matrix) <- paste(samplenames,".2",sep="")
}else{
   rownames(Code_matrix) <- paste(individuals,".2",sep="")
}

HAPMATRIX2 <- Code_matrix

### END OF SECOND HAPLOTYPE MATRIX ################################################################

### THIRD HAPLOTYPE MATRIX ########################################################################
ind_MATRIX <- matrix(,dim(ind)[2],dim(ind)[1])
# fill individual matrix
op <- options(warn = (-1)) # keine Warnmeldungen !
for( xx in 1:dim(ind)[2] ){
  ind_MATRIX[xx,] <- as.integer(substr(ind[,xx], 5, 5)) # 1/1
}
options(op) # Warnmeldungen wieder an !

ind_MATRIX <- ind_MATRIX + 1 

Code_matrix <- matrix(,dim(ind_MATRIX)[1], length(pos)) 
for(xx in 1:dim(MAP_MATRIX)[1]){
 Code_matrix[,xx] <- MAP_MATRIX[xx,][ind_MATRIX[,xx]]
}

Code_matrix <- .Call("code_nucs",Code_matrix)   

gc()
# ----------------------------------------------

# Set the rownames
if(samplenames[1]!=FALSE){
   rownames(Code_matrix) <- paste(samplenames,".3",sep="")
}else{
   rownames(Code_matrix) <- paste(individuals,".3",sep="")
}

HAPMATRIX3 <- Code_matrix

### END OF THIRD HAPLOTYPE MATRIX ################################################################

Code_matrix            <- rbind(HAPMATRIX1,HAPMATRIX2,HAPMATRIX3)

rm(HAPMATRIX1)
rm(HAPMATRIX2)
rm(HAPMATRIX3)
gc()

# modify row Order
NN           <- rownames(Code_matrix)
names(NN)    <- 1:length(NN)
NN           <- sort(NN)
ids          <- as.integer(names(NN))
Code_matrix  <- Code_matrix[ids,,drop=FALSE]
names(NN)    <- NULL
rownames(Code_matrix) <- NN
# ---------------

o_b_j_sub    <- list(matrix=Code_matrix,reference=NaN,positions=pos)

return(o_b_j_sub)
}

# Outsource data # ---------------------------------------
#  dir.create("SNPRObjects")
 
#  o_b_j_sub    <- list(matrix=Code_matrix,reference=NaN,positions=pos)
#  save(o_b_j_sub,file=file.path(getwd(),"SNPRObjects",paste("chr",sep="")))

#---------------------------------------------------------

#  path     <- file.path(getwd(),"SNPRObjects")
#  res      <- readData(path,format="RData",big.data=TRUE,SNP.DATA=TRUE)

#unlink("SNPRObjects", recursive=TRUE)
#return(res)

#}
###############END OF TETRAPLOID DATA




