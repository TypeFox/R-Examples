### GetLinkage

setGeneric("get.recomb", function(object) standardGeneric("get.recomb"))
 setMethod("get.recomb", "GENOME",
 
 function(object){

if(!object@Pop_Recomb$calculated){stop("Statistics have to be calculated first !")}

 res  <- matrix(,object@genelength,1)

  pops <- vector("list",length(object@Pop_Recomb$Populations))
  
  for(xx in 1:length(object@Pop_Recomb$Populations)){
     res[,1]       <- object@RM[,xx]
     
     colnames(res) <- c("Hudson.Kaplan.RM")
     rownames(res) <- object@region.names
     pops[[xx]] <- res 
  }
 
 pops <- as.matrix(pops)
 rownames(pops)  <- paste("pop",1:length(object@Pop_Recomb$Populations))
 colnames(pops)  <- "Recombination"
 
 return(pops)
 
 
 return(res)
 }) 
 
