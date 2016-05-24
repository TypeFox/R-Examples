setGeneric("count.unknowns", function(object) standardGeneric("count.unknowns"))

setMethod("count.unknowns","GENOME",function(object){
  
  region.names                 <- object@region.names
  n.region.names               <- length(region.names)

  if(object@big.data){region.names <- NULL} # because of memory space
 
   
    npops                       <- length(object@populations)     # N Populations
   
 
 #########################################
 # INIT
 #########################################

  nam    <- paste("pop",1:npops)
  init1  <- matrix(0,n.region.names,npops)
  
  missing.nucleotides      <- init1
  missing.frequencies      <- vector("list",n.region.names) # region stats
 
  change    <- object@region.stats
  
#--------------------------------------------------

  # Names ----------------------------------------
  rownames(missing.nucleotides)       <- region.names
  colnames(missing.nucleotides)       <- nam
  # ----------------------------------------------


## PROGRESS #########################
 progr <- progressBar()
#####################################



   
for(xx in 1:n.region.names){

### if Subsites ----------------------------------

bial <- popGetBial(object,xx)

 if(length(bial)!=0){ # if a biallelic position exists
     
   
        
    populations <- object@region.data@populations[[xx]] # Pop 
    

   
    if(length(object@region.data@popmissing[[xx]])!=0){popmissing <- object@region.data@popmissing[[xx]];respop <- (1:npops)[-popmissing]}else{respop <- 1:npops}

   
      res                    <- calc_miss(bial,populations)
    
   # fill detailed Slots --------------------------------#
     missing.frequencies[[xx]]              <- res$miss.freq  
   # ----------------------------------------------------# 
  

    missing.nucleotides[xx,respop]   <- res$miss.nuc
     

  # PROGRESS #######################################################
    progr <- progressBar(xx,n.region.names, progr)
  ###################################################################

 }
}

 change@missing.freqs            <- missing.frequencies
 object@region.stats             <- change
 rm(change)
 gc()
 object@missing.freqs            <- missing.nucleotides

 return(object)
 })


