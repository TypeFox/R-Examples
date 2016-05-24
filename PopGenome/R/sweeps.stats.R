# ------------------------------------------------------------
# Selective Sweeps Statistics
# ------------------------------------------------------------

setGeneric("sweeps.stats", function(object,new.populations=FALSE,subsites=FALSE, freq.table = FALSE, FST = FALSE) standardGeneric("sweeps.stats"))
 setMethod("sweeps.stats", "GENOME",
 function(object,new.populations,subsites, freq.table, FST){

 
 region.names                       <- object@region.names
 n.region.names                     <- length(region.names)
 
# object@Pop_Recomb$sites        <- "ALL"
# object@Pop_Recomb$calculated   <- TRUE
 
 if(missing(new.populations)){
 npops                           <- length(object@populations)
 object@Pop_Sweeps$Populations   <- object@populations
 }else{
 npops                           <- length(new.populations)
 object@Pop_Sweeps$Populations   <- new.populations
 }
 
# bial        <- object@biallelics
# Init
 if(!FST){
  init        <- matrix(NaN,n.region.names,npops)
  CL          <- init
  CLR         <- init
  # Names
  nam                  <- paste("pop",1:npops)
  rownames(CL)         <- region.names
  colnames(CL)         <- nam
  rownames(CLR)        <- region.names
  colnames(CLR)        <- nam
 }else{
  init        <- matrix(NaN,n.region.names,1)
  CL          <- init
  CLR         <- init
 }

 if(!missing(new.populations)){
   NEWPOP      <- TRUE
   populations <- vector("list",npops)
  }else{
   NEWPOP <- FALSE
  } 

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
  # object@Pop_Recomb$sites <- "transitions"   
}

if(subsites=="transversions" & length(bial!=0)){
   transv       <- which(object@region.data@transitions[[xx]]==FALSE)
   bial         <- bial[,transv,drop=FALSE]
  # object@Pop_Recomb$sites <- "transversions"   
} 

if(subsites=="syn" & length(bial!=0)){
   syn        <- which(object@region.data@synonymous[[xx]]==TRUE)
   bial       <- bial[,syn,drop=FALSE]
  # object@Pop_Recomb$sites <- "synonymous"
}
if(subsites=="nonsyn" & length(bial!=0)){
   nonsyn        <- which(object@region.data@synonymous[[xx]]==FALSE)
   bial          <- bial[,nonsyn,drop=FALSE]
  # object@Pop_Recomb$sites <- "nonsynonymous"
}


if(subsites=="intron" & length(bial!=0)){
   intron        <- which(object@region.data@IntronSNPS[[xx]]==TRUE)
   #if(length(intron)==0){
   #       intron <- object@region.data@GeneSNPS[[xx]] & !object@region.data@ExonSNPS[[xx]]	  
   #}
   bial          <- bial[,intron,drop=FALSE]
  # object@Pop_Recomb$sites <- "introns"
}

if(subsites=="utr" & length(bial!=0)){
   utr           <- which(object@region.data@UTRSNPS[[xx]]==TRUE)
   bial          <- bial[,utr,drop=FALSE]
  # object@Pop_Recomb$sites <- "utr"
}

if(subsites=="exon" & length(bial!=0)){
   exon           <- which(object@region.data@ExonSNPS[[xx]]==TRUE)
   bial           <- bial[,exon,drop=FALSE]
  # object@Pop_Recomb$sites <- "exon"
}

if(subsites=="coding" & length(bial!=0)){
   #coding           <- which(!is.na(object@region.data@synonymous[[xx]]))
   coding           <- which(object@region.data@CodingSNPS[[xx]]==TRUE)
   bial             <- bial[,coding,drop=FALSE]
  # object@Pop_Recomb$sites <- "coding"
}

if(subsites=="gene" & length(bial!=0)){
   gene             <- which(object@region.data@GeneSNPS[[xx]])
   bial             <- bial[,gene,drop=FALSE]
  # object@Pop_Recomb$sites <- "gene"
}
}# End if subsites
############### ---------------------------------



  if(length(bial)!=0){ # if a biallelic position exists  
       
    if(NEWPOP){ # wenn neu Populationen definiert
       for(yy in 1:npops){
           if(is.character(new.populations[[yy]])){
              #populations[[yy]] <- match(new.populations[[yy]],rownames(object@DATA[[xx]]@matrix_pol))
              populations[[yy]] <- match(new.populations[[yy]],rownames(bial))
              naids             <- which(!is.na(populations[[yy]]))
              populations[[yy]] <- populations[[yy]][naids]
           }else{
              populations[[yy]] <- new.populations[[yy]]
              ids               <- which(populations[[yy]]>dim(bial)[1])
              if(length(ids)>0){ populations[[yy]] <- populations[[yy]][-ids]}   
           }   
       }
       
       #----------------------#
       temp         <- delNULLpop(populations)
       populations  <- temp$Populations
       popmissing   <- temp$popmissing
       #----------------------#   
       if(length(populations)==0){next} # Keine Population vorhanden
       
    }else{
     populations <- object@region.data@populations[[xx]] # if there is no new population
    }

    
    if(NEWPOP)  {if(length(popmissing)!=0){respop <- (1:npops)[-popmissing]}else{respop <- 1:npops}} # nur die Populationen, die existieren
    if(!NEWPOP) {if(length(object@region.data@popmissing[[xx]])!=0){popmissing <- object@region.data@popmissing[[xx]];respop <- (1:npops)[-popmissing]}else{respop <- 1:npops}}
    

     # important for Coalescent Simulation
     # change@Pop_Linkage[[xx]] <- list(Populations=populations,Outgroup=NULL)
     # ------------------- fill detail slots

     if(!FST){	
     res                       <- complike(bial,populations,freq.table)
     CL[xx,respop]             <- res$CL
    #CLmax[xx,respop]          <- res$CLmax 
     CLR[xx,respop]            <- res$CLR
     # -------------------
     }else{
     res                       <- complike_FST(bial,populations,freq.table)
     CL[xx,] 	               <- res$CL
     CLR[xx,]	               <- res$CLR
     }    

  # PROGRESS #######################################################
    progr <- progressBar(xx,n.region.names, progr)
  ###################################################################

 }
}
 
object@CL    <- CL
object@CLR   <- CLR
#object@CLmax <- CLmax

  return(object)


 })
 





