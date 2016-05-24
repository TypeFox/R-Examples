readSNP <- function(folder, populations=FALSE,outgroup=FALSE, gffpath=FALSE, CHR=FALSE, ref.chr=FALSE, snp.window.size=FALSE, parallized=FALSE, ffpackagebool=TRUE, include.unknown=FALSE)

{


if(!file.exists(folder)){stop("Cannot find folder !")}

if(ref.chr[1]!=FALSE){
 if(!file.exists(ref.chr)){
  stop("Cannot read reference chromosome")
 }
}

if(ref.chr[1]!=FALSE && gffpath==FALSE){
stop("Do not find the GFF/GTF file !")
}


if(file.exists("SNPRObjects")){unlink("SNPRObjects",recursive=TRUE)}
if(file.exists("GFFRObjects")){unlink("GFFRObjects",recursive=TRUE)}


### DEFAULT VALUES for readData()
# include.unknown     <- FALSE
# parallized          <- FALSE
progress_bar_switch <- TRUE
FAST                <- TRUE
big.data            <- TRUE
SNP.DATA            <- TRUE
#### ----------------------------

liste      <- list.files(folder,full.names=TRUE)

# Init 

f.data   <- vector("list",length(liste)) 
f.data2  <- vector("list",length(liste))

if(CHR[1]!=FALSE){
 cat("Read chromosome ", CHR[1], " \n ")
}else{
 cat("Read all chromosomes "," \n ")
}
cat("\n")

individuals     <- character(length(liste))
col.specify     <- c(rep("character",2),"integer",rep("character",2),rep("NULL",4))
col.specify2    <- c(rep("factor",2),"integer",rep("factor",2),rep("NULL",4)) # for read.table.ffdf FIXME

cat("Reading files ... \n")
## PROGRESS #########################
# progr <- progressBar()
#####################################
for (xx in 1:length(liste)){

cat("Reading file: ",liste[[xx]], ": ",xx ," of ", length(liste), "\n" )

 if(CHR[1]!=FALSE){

   CHR             <- as.character(CHR)
   my_pos          <- .Call("find_lines_SNP",liste[[xx]],CHR)
   f.data[[xx]]    <- read.table(file=liste[[xx]],
		      skip=my_pos[1] - 1, nrows=my_pos[2] - my_pos[1],
                      sep="\t",colClasses=col.specify)
 }else{   

   #f.data[[xx]]    <- read.table.ffdf(file=liste[[xx]])
    f.data[[xx]]    <- read.table(file=liste[[xx]],
                       sep="\t",colClasses=col.specify)
 }

 # Split by chromosomes
 chr             <- f.data[[xx]][,2]
 individuals[xx] <- as.character(f.data[[xx]][1,1])

 # chr ids like chr1 or 1
 if(nchar(chr[1])>1){
 chr.pos         <- match(1:5,substr(chr,4,4))
 }else{
 chr.pos         <- match(1:5,chr)
 }

 chromosomes     <- which(!is.na(chr.pos))    # which chromosomes are in the file ?
 chr.pos         <- chr.pos[!is.na(chr.pos)]
 chr.pos         <- c(chr.pos,length(chr))
 
 f.split.data    <- vector("list",length(chr.pos)-1)
 ids             <- 1

 

 for (yy in 1:(length(chr.pos)-1)){
 
   block              <- ids:chr.pos[yy + 1]
   ref_alt            <- as.matrix(f.data[[xx]][block,4:5])
   ref_alt            <- .Call("code_nucs",ref_alt)
   reference.pos      <- f.data[[xx]][block,3]

     #if(ffpackagebool){
     f.split.data[[yy]] <- ff(cbind(ref_alt,reference.pos),dim=c(dim(ref_alt)[1],3)) 
     #}else{
     #f.split.data[[yy]] <- cbind(ref_alt,reference.pos)   
     #}

   ids                <- chr.pos[yy+1]+1

 }

 
f.data2[[xx]] <- f.split.data # seperated by chromosomes

if(CHR[1]==FALSE){
  # delete(f.data[[xx]])
    f.data[[xx]] <- NA
}
# PROGRESS #######################################################
#    progr <- progressBar(xx,length(liste), progr)
###################################################################
}

# RAM freischaufeln !
rm(ref_alt)
rm(reference.pos)
rm(f.data)
rm(f.split.data)
gc()
gc()

# cat("Create Matrix ...\n")
POS          <- NULL

SNP.MATRICES <- vector("list",length(chr.pos)-1)
SNP.POS      <- vector("list",length(chr.pos)-1)

cat("\n")
cat("Prepare ... \n")


 # Create SNP Matrix for each Chromosome
 ## PROGRESS #########################
 progr <- progressBar()
 #####################################

## FOR LOOP START over Chromosomes
for (yy in 1:(length(chr.pos)-1)){
   
   # get positions ------------------------------------
   for(xx in 1:length(liste)){   
       POS              <- c(POS,f.data2[[xx]][[yy]][,3])
       # for saving memory space
       if(!xx%%50){
        POS             <- unique(POS)
       }
   }
   # --------------------------------------------------
  
    
      POS           <- sort(unique(POS))
     
      if(ffpackagebool){     
      SNP.MATRIX    <- ff (NA,dim=c(length(liste),length(POS)), vmode="integer")
      }else{
      SNP.MATRIX    <- matrix(NA,length(liste),length(POS))  
      }

      SNP.POS[[yy]] <- POS
  
   # check what Positions are currently unvisited
   # fill reference first !!!
 
     visited <- vector(,length(POS))
   
     for(xx in 1:length(liste)){ 
    
        # -- neu
	# IDS                <- .Call("my_match_C",f.data2[[xx]][[yy]][,3],POS) 
        # not.na             <- IDS!=-1
        # Keep               <- IDS
        # IDS                <- IDS[not.na]
        # -- neu ende

        # ---alt

         IDS                <- match(f.data2[[xx]][[yy]][,3],POS) 
               
         not.na             <- !is.na(IDS)
         Keep               <- IDS
         IDS                <- IDS[not.na]

        # --- alt ende	

        # welche POS noch nicht besucht 
        IDS                <- IDS[!visited[IDS]]

        if(length(IDS)==0){next}

        IDS2               <- match(IDS,Keep)          
        fill.mat           <- f.data2[[xx]][[yy]][,1][IDS2]
        fill.mat           <- matrix(rep(fill.mat,length(liste)),length(liste),length(fill.mat),byrow=TRUE)        
        SNP.MATRIX[,IDS]   <- fill.mat 
        visited[IDS]       <- TRUE
        if(all(visited)){break}

     }

   rm(fill.mat)
   rm(visited)
   gc()

   # fill subsitution !

    for(xx in 1:length(liste)){ 
 
      # -- neu
      #  IDS                <- .Call("my_match_C",f.data2[[xx]][[yy]][,3],POS) 
      #  not.na             <- IDS!=-1
      #  IDS                <- IDS[not.na]
      # -- neu ende

     # alt ----
       IDS                <- match(f.data2[[xx]][[yy]][,3],POS)
       not.na             <- !is.na(IDS)
       IDS                <- IDS[not.na]
     # alt ende   
             
       SNP.MATRIX[xx,IDS] <- f.data2[[xx]][[yy]][,2] # fill substitution 
     
    }
 
rm(POS)
rm(IDS)
gc()  
  
rownames(SNP.MATRIX)  <- individuals
SNP.MATRICES[[yy]]    <- SNP.MATRIX
rm(SNP.MATRIX)

POS                   <- NULL   

# PROGRESS #######################################################
    progr <- progressBar(yy,(length(chr.pos)-1), progr)
###################################################################

}
### FOR LOOP END over all Chromosomes


#----------------------------------------- Create R Objects ----------------------------------------
 
dir.create("SNPRObjects")

too.big <- FALSE

for ( qq in 1:length(SNP.MATRICES) ) {

 # If there are more than 100000 SNPs

 if(snp.window.size[1]!=FALSE){
 
  too.big     <- TRUE  
  snp.regions <- seq(1,length(SNP.POS[[qq]]),by=snp.window.size)
  snp.regions <- c(snp.regions,length(SNP.POS[[qq]])+1)
  

  for( xx in 1:(length(snp.regions)-1) ){

   CHUNK.MATRIX <- SNP.MATRICES[[qq]][,snp.regions[xx]:(snp.regions[xx+1]-1)]

   o_b_j_sub    <- list(matrix=CHUNK.MATRIX,
                   reference=NaN,positions=SNP.POS[[qq]][snp.regions[xx]:(snp.regions[xx+1]-1)])
   save(o_b_j_sub,file=file.path(getwd(),"SNPRObjects",paste(xx,"chunk_","chr",chromosomes[qq],sep="")))
 
  }  

  # delete.ff(SNP.MATRICES[[qq]])
    SNP.MATRICES[[qq]] <- NA   
    gc()    

 }else{

  o_b_j_sub    <- list(matrix=SNP.MATRICES[[qq]][,],reference=NaN,positions=SNP.POS[[qq]])
  save(o_b_j_sub,file=file.path(getwd(),"SNPRObjects",paste("chr",chromosomes[qq],sep="")))
  # delete.ff(SNP.MATRICES[[qq]])
  SNP.MATRICES[[qq]]  <- NA
  gc()

 }
}

#####--------------- corresponding GFF- File ------------------------#
rm(o_b_j_sub)
gc()


if(gffpath[1]!=FALSE){

cat("\n")
cat("GFF information ...")
cat("\n")
#skip=my_pos[1] - 1, nrows=my_pos[2] - my_pos[1]

 if(CHR[1]!=FALSE){

   my_pos      <- .Call("find_lines_GFF",gffpath,CHR)
   gff.table   <- read.table(gffpath,sep="\t",colClasses=c(rep("character",3),rep("integer",2),rep("character",2),"character","NULL"),
                  skip = my_pos[1] - 1, nrows = my_pos[2] - my_pos[1] + 1)
 }else{

   gff.table   <- read.table(gffpath,sep="\t",colClasses=c(rep("character",3),rep("integer",2),rep("character",2),"character","NULL"))

 }

 dir.create("GFFRObjects")

### sort GFF file
POS        <- gff.table[,4]
names(POS) <- 1:length(POS)
POS        <- sort(POS)
ids        <- as.numeric(names(POS))
gff.table  <- gff.table[ids,,drop=FALSE]
## vorletzte Zeile integer draus wegen GFF3 "." -> 0
vorl       <- gff.table[,8]
punkte     <- which(vorl==".")
if(length(punkte)>0){
  vorl[punkte]  <- 0
  gff.table[,8] <- as.integer(vorl) 
}else{
  gff.table[,8] <- as.integer(gff.table[,8])
}

#-------------------------
seq.names  <- gff.table[,1]


# If the data for one chromosome is too big !!!
if(too.big){
 
 # - works only for one chromosome
 # ---------------------------------

 IDSS       <- grep(CHR,seq.names)
 gff.table  <- gff.table[IDSS,,drop=FALSE]

 # for every chunk of ONE chromosome
 for( xx in 1:(length(snp.regions)-1) ){
    # get the region 
    pos1         <- SNP.POS[[1]][snp.regions[xx]]
    pos2         <- SNP.POS[[1]][(snp.regions[xx+1]-1)]
    #cat(pos1,"-",pos2,"  \n")
    SUB_GFF      <- split.GFF(gff.table,pos1,pos2)
    o_b_j_sub    <- SUB_GFF
    save(o_b_j_sub,file=file.path(getwd(),"GFFRObjects",paste(xx,"chunk_","chr",CHR[1],sep="")))  
 }

}else{
 
 for(i in chromosomes){
    pos          <- grep(i,seq.names)
    o_b_j_sub    <- gff.table[pos,]
    save(o_b_j_sub,file=file.path(getwd(),"GFFRObjects",paste("chr",i,sep="")))
 } 
}

gffpath  <- file.path(getwd(),"GFFRObjects")
   
}

######### --------------- End of get GFF ---------------------------------#

path     <- file.path(getwd(),"SNPRObjects")
cat("\n")
cat("Calculation ... \n")
res      <- readData(path,populations=populations,outgroup=outgroup,include.unknown=include.unknown,gffpath=gffpath,format="RData",
            parallized=parallized,progress_bar_switch=progress_bar_switch,
            FAST=FAST,big.data=big.data,SNP.DATA=SNP.DATA)

if(too.big){
 res              <- concatenate_to_whole_genome(res,length(res@region.names))
 res@region.names <- paste("chr",CHR[1],sep="")
}


# Loeschen des Verzeichnisses

unlink("SNPRObjects", recursive=TRUE)
if(gffpath[1]!=FALSE){
unlink("GFFRObjects", recursive=TRUE)
}


if(ref.chr[1]!=FALSE){
cat("\n")
cat("Synonymous positions ...")
cat("\n")
res <- set.synnonsyn(res,ref.chr=ref.chr)
}

gc()

return(res)

}


