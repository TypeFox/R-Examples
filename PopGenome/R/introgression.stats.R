setGeneric("introgression.stats", function(object,subsites=FALSE,do.D=TRUE) standardGeneric("introgression.stats"))

setMethod("introgression.stats","GENOME",function(object,subsites,do.D){
  
  region.names                 <- object@region.names
  n.region.names               <- length(region.names)
  if(object@big.data){region.names <- NULL} # because of memory space
 
    NEWPOP <- FALSE
    npops                       <- length(object@populations)     # total number of populations
     
   
# init
  D      <- matrix(NaN,n.region.names,1)
  f      <- matrix(NaN,n.region.names,1)
#--------------------------------------------------

  # Names ----------------------------------------
  rownames(D)       <- region.names
  colnames(D)       <- "D statistic"
  rownames(f)       <- region.names
  colnames(f)       <- "f statistic"
  # ----------------------------------------------


## PROGRESS #########################
 progr <- progressBar()
#####################################
   
for(xx in 1:n.region.names){

### if Subsites ----------------------------------

bial <- popGetBial(object,xx)

if(subsites[1]!=FALSE){

if(subsites=="transitions" & length(bial!=0)){
   tran       <- which(object@region.data@transitions[[xx]]==TRUE)
   bial       <- bial[,tran,drop=FALSE]
  # object@Pop_FSTH$sites <- "transitions"   
}

if(subsites=="transversions" & length(bial!=0)){
   transv       <- which(object@region.data@transitions[[xx]]==FALSE)
   bial         <- bial[,transv,drop=FALSE]
  # object@Pop_FSTH$sites <- "transversions"   
} 

if(subsites=="syn" & length(bial!=0)){
   syn        <- which(object@region.data@synonymous[[xx]]==TRUE)
   bial       <- bial[,syn,drop=FALSE]
  # object@Pop_FSTH$sites <- "synonymous"
}
if(subsites=="nonsyn" & length(bial!=0)){
   nonsyn        <- which(object@region.data@synonymous[[xx]]==FALSE)
   bial          <- bial[,nonsyn,drop=FALSE]
  # object@Pop_FSTH$sites <- "nonsynonymous"
}

if(subsites=="intron" & length(bial!=0)){
   intron        <- which(object@region.data@IntronSNPS[[xx]]==TRUE)
   #if(length(intron)==0){
   #       intron <- object@region.data@GeneSNPS[[xx]] & !object@region.data@ExonSNPS[[xx]]	  
   #}
   bial          <- bial[,intron,drop=FALSE]
  # object@Pop_Linkage$sites <- "introns"
}

if(subsites=="utr" & length(bial!=0)){
   utr           <- which(object@region.data@UTRSNPS[[xx]]==TRUE)
   bial          <- bial[,utr,drop=FALSE]
  # object@Pop_FSTH$sites <- "utr"
}

if(subsites=="exon" & length(bial!=0)){
   exon           <- which(object@region.data@ExonSNPS[[xx]]==TRUE)
   bial           <- bial[,exon,drop=FALSE]
  # object@Pop_FSTH$sites <- "exon"
}

if(subsites=="coding" & length(bial!=0)){
   #coding           <- which(!is.na(object@region.data@synonymous[[xx]])==TRUE)
   coding           <- which(object@region.data@CodingSNPS[[xx]]==TRUE)
   bial             <- bial[,coding,drop=FALSE]
  # object@Pop_FSTH$sites <- "coding"
}

if(subsites=="gene" & length(bial!=0)){
   gene             <- which(object@region.data@GeneSNPS[[xx]]==TRUE)
   bial             <- bial[,gene,drop=FALSE]
  # object@Pop_FSTH$sites <- "gene"
}

if(subsites=="intergenic"){
  intron        <- which(object@region.data@IntronSNPS[[xx]]==TRUE)
   if(length(intron)==0){
     intron <- !object@region.data@ExonSNPS[[xx]]	  
   }

  utr            <- object@region.data@UTRSNPS[[xx]]
  exon           <- object@region.data@ExonSNPS[[xx]]
  gene           <- object@region.data@GeneSNPS[[xx]]
  coding         <- !is.na(object@region.data@synonymous[[xx]])  

  inter          <- !(intron|utr|exon|gene|coding)
  bial           <- bial[,inter,drop=FALSE]
  #object@Pop_FSTH$sites <- "intergenic"
}
}# End if subsites
############### ---------------------------------


 if(length(bial)!=0){ # if a biallelic position exists
           
    populations <- object@region.data@populations[[xx]] 
    outgroup    <- object@region.data@outgroup[[xx]]  

    if(do.D){
    res       <- calc_D(bial, populations, outgroup) 
    D[xx]     <- res$D
    f[xx]     <- res$f

    }

  ## PROGRESS #######################################################
    progr <- progressBar(xx,n.region.names, progr)
  ###################################################################

 }
}

 object@D  <- D
 object@f  <- f
 
 return(object)
 })


