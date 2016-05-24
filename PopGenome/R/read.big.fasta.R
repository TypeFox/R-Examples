read.big.fasta <- function(filename,populations=FALSE,outgroup=FALSE,window=2000,SNP.DATA=FALSE,include.unknown=FALSE,
parallized=FALSE,FAST=FALSE,big.data=TRUE){


if(!file.exists(filename)){stop("Cannot find file !")}

if(file.exists("FASTARObjects")){unlink("FASTARObjects",recursive=TRUE)}


### DEFAULT VALUES for readData()
progress_bar_switch <- TRUE
gffpath=FALSE
#### ----------------------------


info   <- .Call("get_dim_fasta",filename)

if(length(info)==0){
stop("Can not open file \n")
}

dim.info   <- info[[1]]
ind.names  <- info[[2]]

#print(dim.info)
#print(ind.names)
#Create ff matr
if(big.data){
MAT  <- ff(5,dim=c(dim.info[1],dim.info[2]))
}
if(!big.data){
MAT  <- matrix(5,nrow=dim.info[1],ncol=dim.info[2])
}
#MAT   <- matrix(0,dim.info[1],dim.info[2])

 fasta_handle <- .Call("FASTA_open",filename,"r")

  for (xx in 1:length(ind.names) ){
     #print(xx)
     #print(dim.info[2])
     #vec <- .Call("get_ind_fasta",filename,xx,dim.info[2])    
     #print(length(vec))
     vec      <- .Call("FASTA_getNextIndividual", fasta_handle, dim.info[2])
     MAT[xx,] <- vec
     cat(xx," of ", dim.info[1], " individuals","\n")
     
  }

rm(fasta_handle)
gc()

rownames(MAT) <- ind.names

# Outsource Chunks
chunks <- seq(1,dim.info[2],by=window)
if(chunks[length(chunks)]!=dim.info[2]){chunks <- c(chunks,dim.info[2])}

 #if(.Platform$OS.type=="windows"){
 # shell("md FASTARObjects")
 #}else{
 # system("mkdir FASTARObjects")
 #}

 dir.create("FASTARObjects")

## Outsource first matrix
chunk.matrix <- MAT[,chunks[1]:chunks[2]]
o_b_j_sub    <- list(matrix=chunk.matrix,reference=NaN,positions=NaN)
save(o_b_j_sub,file=file.path(getwd(),"FASTARObjects",paste("chunk",1,sep="")))

for (xx in 2:(length(chunks))-1){

 chunk.matrix <- MAT[,(chunks[xx]+1):chunks[xx+1]]
 o_b_j_sub    <- list(matrix=chunk.matrix,reference=NaN,positions=NaN)
 save(o_b_j_sub,file=file.path(getwd(),"FASTARObjects",paste("chunk",xx,sep="")))

}

path     <- file.path(getwd(),"FASTARObjects")

res      <- readData(path,populations=populations,outgroup=outgroup,include.unknown=include.unknown,gffpath=gffpath,format="RData",
            parallized=parallized,progress_bar_switch=progress_bar_switch,
            FAST=FAST,big.data=big.data,SNP.DATA=SNP.DATA)

# Loeschen des Verzeichnisses


unlink("FASTARObjects", recursive=TRUE)


#if(.Platform$OS.type=="windows"){
#  shell(paste("rmdir ",file.path(getwd(),"FASTARObjects")," /s /q",sep=" "))
#  if(gffpath[1]!=FALSE){
#  shell(paste("rmdir ",file.path(getwd(),"GFFRObjects")," /s /q",sep=" "))
#  }  
#}else{
#  system(paste("rm -r ",file.path(getwd(),"FASTARObjects"),sep=" "))
#  if(gffpath[1]!=FALSE){
#  system(paste("rm -r ",file.path(getwd(),"GFFRObjects"),sep=" "))
#  }
#}


return(res)

}
