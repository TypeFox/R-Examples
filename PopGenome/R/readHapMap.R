readHapMap <- function(folder, hap_gffpath, populations=FALSE,outgroup=FALSE){


# Default values
SNP.DATA            <- TRUE
big.data            <- TRUE
FAST                <- TRUE
progress_bar_switch <- TRUE
parallized          <- FALSE
gffpath             <- FALSE
include.unknown     <- FALSE
chromosomes         <- 1:22 # FIXME
#---------------------


liste      <- list.files(folder,full.names=TRUE)
liste2     <- list.files(folder)

gff_liste  <- list.files(hap_gffpath,full.names=TRUE)
gff_liste2 <- list.files(hap_gffpath)


rbind_list  <- vector("list",length(chromosomes)) # Hapmap
rbind_list2 <- vector("list",length(chromosomes)) # gff

 yy <- 1
 for(xx in chromosomes){

 chr                <- paste("chr",xx,sep="") 
 what               <- grep(chr,liste2,fixed=TRUE)
 what_gff           <- grep(chr,gff_liste2,fixed=TRUE)
 rbind_list[[yy]]   <- liste[what] 
 rbind_list2[[yy]]  <- gff_liste[what_gff]
 yy <- yy + 1

}

for(xx in 1:length(chromosomes)){

# 1e+08 ist zu gross ca.
o_b_j_       <- create_SNP_matrix(rbind_list[[xx]],rbind_list2[[xx]])

# hier koennte man dann anfangen zu splitten
rows <- dim(o_b_j_$matrix)[1]
cols <- dim(o_b_j_$matrix)[2]

if ((rows*cols)>1e+08){

chunk        <- 10000
n.blocks     <- ceiling(cols/10000)
blocks       <- split(1:cols,1:n.blocks)    

for(yy in 1:n.blocks){

o_b_j_sub    <- list(matrix=o_b_j_$matrix[,blocks[[yy]]],reference=NaN,positions=o_b_j_$positions[blocks[[yy]]])

 #if(.Platform$OS.type=="windows"){
 # shell("md HapMapRObjects")
 #}else{
 # system("mkdir HapMapRObjects")
 #}
 
 dir.create("HapMapRObjects")


  save(o_b_j_sub,file=file.path(getwd(),"HapMapRObjects",paste("chr",chromosomes[xx],"chunk",yy,sep="")))
 

}

}else{ # end of splitting


o_b_j_sub    <- list(matrix=o_b_j_$matrix[],reference=NaN,positions=o_b_j_$positions[])

 #if(.Platform$OS.type=="windows"){
 # shell("md HapMapRObjects")
 #}else{
 # system("mkdir HapMapRObjects")
 #}

 dir.create("HapMapRObjects")
 
  save(o_b_j_sub,file=file.path(getwd(),"HapMapRObjects",paste("chr",chromosomes[xx],sep="")))
 



}
rm(o_b_j_)
gc() #?=

}
path     <- file.path(getwd(),"HapMapRObjects")
res      <- readData(path,populations=populations,outgroup=outgroup,include.unknown=include.unknown,gffpath=gffpath,format="RData",
            parallized=parallized,progress_bar_switch=progress_bar_switch,
            FAST=FAST,big.data=big.data,SNP.DATA=SNP.DATA)

# Loeschen des Verzeichnisses

#if(.Platform$OS.type=="windows"){
#  shell(paste("rmdir ",file.path(getwd(),"HapMapRObjects")," /s /q",sep=" "))
# }else{
#  system(paste("rm -r ",file.path(getwd(),"HapMapRObjects"),sep=" "))
# }

unlink("HapMapRObjects",recursive=TRUE)


return(res)
}

## SUBFUNCTIONS######################################

create_SNP_matrix <- function(files,gff_files){

# for one chromosome

matrices    <- vector("list",length(files))
pos         <- vector("list",length(files))
positions   <- NULL
individuals <- NULL
individuals.names <- NULL
reference   <- vector("list",length(gff_files))

for(xx in 1:length(files)){
               BB  <- parse_HapMap(files[xx])
   matrices[[xx]]  <- BB$matrix
   pos[[xx]]       <- BB$positions
   #print(dim(matrices[[xx]])[2])
   #print(length(pos[[xx]]))
   #print(length(unique(pos[[xx]])))
				
#print("!!")	
#print(pos[[xx]][length(pos[[xx]])])
						
   positions          <- c(positions,BB$positions)
   individuals        <- c(individuals,dim(BB$matrix)[1])
   individuals.names  <- c(individuals.names,rownames(BB$matrix))           
   reference[[xx]]    <- get_ref_alleles(gff_files[xx])[-BB$dupids] 
#   print(length(reference[[xx]]))  
}

positions <- sort(unique(positions))
Mrow      <- sum(individuals)
Mcol      <- length(positions)
####---- FF ------------
MATRIX    <- ff(matrix(NA,Mrow,Mcol),dim=c(Mrow,Mcol),vmode="integer")
########################

#print(dim(MATRIX))
# Fill Matrix

#print("individuals")
#print(individuals)
rows     <- 1:individuals[1]
pop.rows <- vector("list",length(files))
for(xx in 1:length(files)){
 pop.rows[[xx]] <- rows
# print(length(rows))
 ids    <- is.element(positions,pos[[xx]])
# print(length(pos[[xx]]))
# print(dim(matrices[[xx]]))
# print(length(ids))
# print(length(which(ids)))
 MATRIX[rows,ids] <- matrices[[xx]] 
 if(xx<length(files)){
 rows <- (rows[length(rows)]+1):(rows[length(rows)]+individuals[xx+1]) 
 }
}

cat("Fill NA positions with reference allele ...\n")
# Fill NA positions
# na.ids  <- which(is.na(MATRIX),arr.ind=TRUE)
# na.cols <- unique(na.ids[,2])
#---FF----
ja_na   <- function(x){return(any(is.na(x)))}
na.ids  <- ffapply(X=MATRIX,MARGIN=2,AFUN="ja_na",RETURN=TRUE,CFUN="list")
na.cols <- which(unlist(na.ids))
#---FF----

na.pos  <- positions[na.cols]

# use pop.rows
# Flaschenhals.de

## PROGRESS #########################
 progr <- progressBar()
#####################################

for(xx in 1:length(na.pos)){
 for(yy in 1:length(files)){
    id <- match(na.pos[xx],pos[[yy]])
    if(!is.na(id)){
      column  <- MATRIX[,na.cols[xx]]
      rowws   <- which(is.na(column))
      MATRIX[rowws,na.cols[xx]] <- reference[[yy]][id]
    break
    }
 }
# PROGRESS #######################################################
    progr <- progressBar(xx,length(na.pos), progr)
################################################################### 
}

#iii <<- 1
#apply(MATRIX,2,function(x){
#nana <- which(is.na(x))
#if(length(nana)>0){
#  na.pos <- positions[iii]
#  # find reference allele  
#  kkk <<- 1
#  breakpoint  <<- TRUE
#  lapply(reference,function(y){
#    id <- match(na.pos,y)
#    if(!is.na(id) & breakpoint){
#     MATRIX[nana,iii] <<- reference[[kkk]][id]  
#     breakpoint <<- FALSE
#    }
#   
#  kkk <<- kkk + 1  
#  })	
#}
#iii     <<- iii + 1
#counter <-  iii 
# PROGRESS #######################################################
#    progr <- progressBar(counter,dim(MATRIX)[2], progr)
################################################################### 
#})

#### FF -------
dimnames(MATRIX) <- list(individuals.names,NULL) # FFFFFFFFFFFFF
# print(MATRIX[,1:2])
return(list(matrix=MATRIX,positions=positions))
}

get_ref_alleles <- function(XXX){
pos <- read.table(XXX,sep=";",colClasses=c("NULL","NULL","character","NULL","NULL","NULL"))
pos <- as.matrix(pos)
pos <- substr(pos,11,11)
pos <- .Call("code_nucs",pos)
return(pos)
}

