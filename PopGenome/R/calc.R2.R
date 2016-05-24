setGeneric("calc.R2", function(object,subsites=FALSE,lower.bound=0, upper.bound=1) standardGeneric("calc.R2"))
 setMethod("calc.R2", "GENOME",
 function(object,subsites,lower.bound,upper.bound){

 region.names                       <- object@region.names
 n.region.names                     <- length(region.names)
 if(object@big.data){region.names <- NULL} # because of memory space
 
 
 object@Pop_Linkage$sites        <- "ALL"
 object@Pop_Linkage$calculated   <- TRUE
 
 npops                           <- length(object@populations)
 object@Pop_Linkage$Populations  <- object@populations

  
# Init
 linkage.disequilibrium <- vector("list",n.region.names)
 change      <- object@region.stats
 

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
   #object@Pop_Linkage$sites <- "transitions"   
}

if(subsites=="transversions" & length(bial!=0)){
   transv       <- which(object@region.data@transitions[[xx]]==FALSE)
   bial         <- bial[,transv,drop=FALSE]
   #object@Pop_Linkage$sites <- "transversions"   
} 

if(subsites=="syn" & length(bial!=0)){
   syn        <- which(object@region.data@synonymous[[xx]]==TRUE)
   bial       <- bial[,syn,drop=FALSE]
   #object@Pop_Linkage$sites <- "synonymous"
}
if(subsites=="nonsyn" & length(bial!=0)){
   nonsyn        <- which(object@region.data@synonymous[[xx]]==FALSE)
   bial          <- bial[,nonsyn,drop=FALSE]
   #object@Pop_Linkage$sites <- "nonsynonymous"
}


if(subsites=="intron" & length(bial!=0)){
   intron        <- which(object@region.data@IntronSNPS[[xx]]==TRUE)
   #if(length(intron)==0){
   #       intron <- object@region.data@GeneSNPS[[xx]] & !object@region.data@ExonSNPS[[xx]]	  
   #}
   bial          <- bial[,intron,drop=FALSE]
   #object@Pop_Linkage$sites <- "introns"
}

if(subsites=="utr" & length(bial!=0)){
   utr           <- which(object@region.data@UTRSNPS[[xx]]==TRUE)
   bial          <- bial[,utr,drop=FALSE]
   #object@Pop_Linkage$sites <- "utr"
}

if(subsites=="exon" & length(bial!=0)){
   exon           <- which(object@region.data@ExonSNPS[[xx]]==TRUE)
   bial           <- bial[,exon,drop=FALSE]
   #object@Pop_Linkage$sites <- "exon"
}

if(subsites=="coding" & length(bial!=0)){
   #coding           <- which(!is.na(object@region.data@synonymous[[xx]]))
   coding           <- which(object@region.data@CodingSNPS[[xx]]==TRUE)
   bial             <- bial[,coding,drop=FALSE]
   #object@Pop_Linkage$sites <- "coding"
}

if(subsites=="gene" & length(bial!=0)){
   gene             <- which(object@region.data@GeneSNPS[[xx]])
   bial             <- bial[,gene,drop=FALSE]
   #object@Pop_Linkage$sites <- "gene"
}
}# End if subsites
############### ---------------------------------



  if(length(bial)!=0){ # if a biallelic position exists  
       
     populations <- object@region.data@populations[[xx]] # if there is no new population
   
       if((lower.bound!=0) | (upper.bound !=1)){
        # Get the frequencies
        freq <- jointfreqdist(bial,list(unlist(populations))) 
	freq <- freq$jfd        	
        sub       <- (freq >= lower.bound) & (freq <=upper.bound)
        bial      <- bial[,sub,drop=FALSE]
       }

     res                            <- intern.calc.R2(bial,populations)
     linkage.disequilibrium[[xx]]   <- res$res
     
  # PROGRESS #######################################################
    progr <- progressBar(xx,n.region.names, progr)
  ###################################################################

 }
}

 change@linkage.disequilibrium <- linkage.disequilibrium
 object@region.stats <- change 
  
return(object)
 })
 
