intern.calc.R2 <- function(matrix_pol,populations){

### NO NA in DATA

## NOTE: Gaps in the biallelic matrix should be deleted
npops        <- length(populations)
sitelength   <- dim(matrix_pol)[2]

if(sitelength<2){
 return(list(res=as.matrix(NaN)))
}

if(length(colnames(matrix_pol))==0){
 colnames(matrix_pol) <- 1:sitelength
}

 segsites <- get_segsites_FAST(matrix_pol,populations) # positions of the segsites of each population


reslist <- vector("list", npops)

for(xx in 1:npops){

 popmatrix        <- matrix_pol[populations[[xx]],segsites[[xx]],drop=FALSE]
 n.segsites.pop   <- length(segsites[[xx]])
 

 if(n.segsites.pop<=1){next}

  bialsites    <- as.numeric(colnames(matrix_pol))
  EINSEN       <- colSums(popmatrix, na.rm=TRUE)
  NULLEN       <- colSums(popmatrix==0, na.rm=TRUE)
  res          <- .Call("R2_C_plus", popmatrix, EINSEN, NULLEN, bialsites)  
 
 R2  <- res[[1]]
 M   <- res[[2]]

 P   <- apply(M,1,function(x){return(fisher.test( matrix(x,ncol=2,byrow=TRUE) )$p.value ) 
             })
 
 d_dist <- abs(res[[3]])
 res <- rbind(R2,P,d_dist)
 rownames(res) <-c("R2","P","Distance")

reslist[[xx]] <- res


}# End of iteration over pops

 reslist           <- as.matrix(reslist)
 rownames(reslist) <- paste("pop",1:npops)
 colnames(reslist) <- "Linkdisequilibrium"
 
 return(list(res=reslist))

}



