setGeneric("getVolumeUD", function(x,...){standardGeneric("getVolumeUD")})
setMethod("getVolumeUD", 
          signature=c(x=".UDStack"),
          definition=function(x,...){
		  stack(lapply(split(x), getVolumeUD))
	  })
setMethod("getVolumeUD", 
          signature=c(x=".UD"),
          definition=function(x,...){
            transf <- function(nr){
              rank <- (1:length(values(nr)))[rank(values(nr))]
              values(nr) <- 1 - cumsum(sort(values(nr)))[rank]
              return(as(nr,'RasterLayer'))
            }            
            if(length(c(x,...))==1) {return(transf(x))} else {return(lapply(list(x, ...), transf))}
            })
