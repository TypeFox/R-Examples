setGeneric("mult.linkage.stats", function(object,lower.bound=0,upper.bound=1,pairs=FALSE) standardGeneric("mult.linkage.stats"))
 setMethod("mult.linkage.stats", "GENOME",

 function(object,lower.bound,upper.bound,pairs){


 # define an evironment
 multGLOBAL <- new.env()

 if(pairs[1]==FALSE){
  pairs      <- combn(length(object@region.names),2)
 }

 n.pairs    <- dim(pairs)[2]

 multGLOBAL$res        <- vector("list",n.pairs)
 multGLOBAL$iter       <- 1

#### NAMES ----------------------------------------
pp <- pairs
nn <- paste("",pp[1,1],"/",pp[2,1],sep="")
if(dim(pp)[2]>1){ # more than 2 sites
 for(yy in 2:dim(pp)[2]){
    m <- paste("",pp[1,yy],"/",pp[2,yy],sep="")
    nn <- c(nn,m)
 }
}#END if
#### ---------


exx <- apply(pairs,2,function(xx){
	

cat("region", xx[1] ,"vs ", xx[2], "\n")
	
	bial1 <- popGetBial(object,xx[1])
	bial2 <- popGetBial(object,xx[2])

	pop1 <- object@region.data@populations[[xx[1]]]
	pop2 <- object@region.data@populations[[xx[2]]]

        if(length(bial1)==0){multGLOBAL$iter <- multGLOBAL$iter + 1;return(0)}
	if(length(bial2)==0){multGLOBAL$iter <- multGLOBAL$iter + 1;return(0)}

	# Get the frequencies
        freq1 <- jointfreqdist(bial1,pop1) 
        freq2 <- jointfreqdist(bial2,pop2) 
 
	freq1 <- freq1$jfd        
	freq2 <- freq2$jfd

        subsites1 <- (freq1 >= lower.bound) & (freq1 <=upper.bound)
	subsites2 <- (freq2 >= lower.bound) & (freq2 <=upper.bound)
       	

        bial1 <- bial1[,subsites1,drop=FALSE]
        bial2 <- bial2[,subsites2,drop=FALSE]

 	multGLOBAL$res[[multGLOBAL$iter]]   <- pair_linkdisequ_FAST(bial1,bial2,pop1,pop2)
            

        multGLOBAL$iter          <- multGLOBAL$iter + 1

 })

res           <- as.matrix(multGLOBAL$res)
rownames(res) <- nn
colnames(res) <- "region-pairwise Linkage"
object@mult.Linkage <- res

return(object)

})
