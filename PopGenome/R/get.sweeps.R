### GetSweeps

setGeneric("get.sweeps", function(object) standardGeneric("get.sweeps"))
 setMethod("get.sweeps", "GENOME",
 
 function(object){

# if(!object@Pop_Recomb$calculated){stop("Statistics have to be calculated first !")}

 res  <- matrix(,object@genelength,2)

  pops <- vector("list",length(object@Pop_Sweeps$Populations))
  
  for(xx in 1:length(object@Pop_Sweeps$Populations)){
     res[,1]       <- object@CL[,xx]
     res[,2]       <- object@CLmax[,xx]
     
     colnames(res) <- c("Nielsen CL","Nielsen CLmax")
     rownames(res) <- object@region.names
     pops[[xx]]    <- res 

  }
 
 pops            <- as.matrix(pops)
 rownames(pops)  <- paste("pop",1:length(object@Pop_Sweeps$Populations))
 colnames(pops)  <- "Selective Sweeps"
 
 return(pops)
 
 
 return(res)
 }) 
 
