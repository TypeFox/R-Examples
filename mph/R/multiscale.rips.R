multiscale.rips.ipca <- function(X, epsGMRA, maxD, single=FALSE){
   n <- ncol(X)
   m <- nrow(X)
   res <- .Call("multiscale_rips_ipca", as.double(t(X)), m, n,
        as.double(epsGMRA), as.integer(maxD), as.integer(single) )  
   
   colnames(res) <- c( "birth", "death", "scale", "dimension" );
   res[,4] = res[,4]-2
   res 
}
