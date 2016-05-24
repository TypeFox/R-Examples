setGeneric("region.as.fasta", function(object, region.id= FALSE, filename = FALSE, type=1, ref.chr=FALSE) standardGeneric("region.as.fasta"))
setMethod("region.as.fasta", "GENOME",
function(object,region.id,filename,type,ref.chr){

if(region.id[1]==FALSE){
 stop("Define a region !")
}


if(type==1){

bial    <- popGetBial(object,region.id)
subst   <- object@region.data@biallelic.substitutions[[region.id]]
minor   <- subst[1,]
mayor   <- subst[2,]

 for(xx in 1:dim(bial)[2]){

  vek       <- bial[,xx]
  ids.minor <- vek==1
  ids.mayor <- vek==0

  bial[ids.minor,xx] <- minor[xx]
  bial[ids.mayor,xx] <- mayor[xx]

 }

 

 number     <- c(1,1,1,1,2,2,3,3,4,4,5,5,5,6)
 nuc        <- c("T","t","U","u","C","c","G","g","A","a","N","n","?","-")


 bial <- apply(bial,1,function(x){return(nuc[match(x,number)])})
 bial <- t(bial)

APE_write.dna(bial,file=filename,colsep="",format="fasta") # This function is from the APE package on CRAN 

return(bial)
}


if(type==2){

#  Init whole MATRIX
## Reading the reference chromosome
file.info  <- .Call("get_dim_fasta",ref.chr)
CHR        <- .Call("get_ind_fasta",ref.chr,1,file.info[[1]][2])
bial.sites <-  object@region.data@biallelic.sites[[region.id]]
s_tart     <-  bial.sites[1]
e_end      <-  bial.sites[length(bial.sites)]     

bial      <-  popGetBial(object,region.id)
FILLnuc   <-  CHR[s_tart:e_end]
# whole Matrix
RETMAT    <-  matrix(rep(FILLnuc,dim(bial)[1]),nrow=dim(bial)[1],ncol=length(FILLnuc),byrow=TRUE)


subst   <- object@region.data@biallelic.substitutions[[region.id]]
minor   <- subst[1,]
mayor   <- subst[2,]


# fill the SNPS
 for(xx in 1:dim(bial)[2]){

  vek       <- bial[,xx]
  ids.minor <- vek==1
  ids.mayor <- vek==0

  bial[ids.minor,xx] <- minor[xx]
  bial[ids.mayor,xx] <- mayor[xx]

 }

 number     <- c(1,1,1,1,2,2,3,3,4,4,5,5,5,6)
 nuc        <- c("T","t","U","u","C","c","G","g","A","a","N","n","?","-")

 bial     <- apply(bial,1,function(x){return(nuc[match(x,number)])})
 bial     <- t(bial)

ind.names <- rownames(bial)

 RETMAT   <- apply(RETMAT,1,function(x){return(nuc[match(x,number)])})
 RETMAT   <- t(RETMAT)

fillids   <- match(bial.sites,s_tart:e_end)

RETMAT[,fillids] <- bial

rownames(RETMAT) <- ind.names

# if a region is specified via splitting.data concatenate left and right-hand nucleotides
if( length(grep("-",object@region.names[region.id]))!=0 ){

	leftrightPOS <- as.numeric(strsplit(object@region.names[region.id]," - ")[[1]])
        leftPOS      <- leftrightPOS[1] 
	rightPOS     <- leftrightPOS[2]
        # s_tart 
	# ---------------------------------------------------------------------------------
        # left hand matrix
	leftRegion   <- leftPOS:(s_tart-1)
        leftNUCS     <- CHR[leftRegion]
        # convert to character
	leftNUCS     <- nuc[match(leftNUCS,number)] 
        leftNUCS     <- rep(leftNUCS,length(ind.names))
        leftNUCS     <- matrix(leftNUCS,nrow=length(ind.names),ncol=length(leftRegion),byrow=TRUE)
        #e_end
	# ----------------------------------------------------------------------------------
	# right hand matrix
	rightRegion  <- (e_end+1):rightPOS
        rightNUCS    <- CHR[rightRegion]
        # convert to character
	rightNUCS    <- nuc[match(rightNUCS,number)] 
        rightNUCS    <- rep(rightNUCS,length(ind.names))
        rightNUCS    <- matrix(rightNUCS,nrow=length(ind.names),ncol=length(rightRegion),byrow=TRUE)

# concatenate matrices

RETMAT <- cbind(leftNUCS,RETMAT,rightNUCS)

}


APE_write.dna(RETMAT,file=filename,colsep="",format="fasta") # This function is from the ape package on CRAN 

return(RETMAT)


}

})# End of function

