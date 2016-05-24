calc_miss <- function(matrix_pol,populations){

# Only one polymorphic site
if(is.vector(matrix_pol)){
   matrix_pol <- as.matrix(matrix_pol)
   warning("#---------------> Only one polymorphic site <-------------------#")
}

npops         <- length(populations)

miss.nuc      <- numeric(npops)
miss.freq     <- matrix(0,npops,dim(matrix_pol)[2])


 for(xx in 1:npops){

 popbial  <- matrix_pol[populations[[xx]],,drop=FALSE]
 All      <- dim(popbial)[1]*dim(popbial)[2]
 n.M      <- sum(is.na(popbial))
 miss.nuc[xx] <- n.M/All

 # Now for each SNP 
 miss.freq <- apply(popbial,2,function(x){
	   
           All <- length(x)
           n.M <- sum(is.na(x))
           return(n.M/All)
 })
 

 }# end for over pops


return(list(miss.nuc=miss.nuc,miss.freq=miss.freq))

}


