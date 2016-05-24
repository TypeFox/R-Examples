pair_linkdisequ_FAST <- function(bial1, bial2, populations1, populations2){


if(length(colnames(bial1))>1){

# only need for Distance !!
combinations    <- NULL
vek             <- 1:dim(bial1)[2]
vek2            <- 1:dim(bial2)[2]

 for(xx in vek){
  val <- rep(xx,dim(bial2)[2])
  mat <- rbind(val,vek2)
  combinations <- cbind(combinations,mat) 
 }

}
# --------------------------------

npops    <- length(populations1)
reslist  <- vector("list",npops)

for(yy in 1:npops){

bial1_sub <- bial1[populations1[[yy]],,drop=FALSE]
bial2_sub <- bial2[populations2[[yy]],,drop=FALSE]

colsums1 <- colSums(bial1_sub)
colsums2 <- colSums(bial2_sub)

res <- .Call("R2_between_C",bial1_sub,colsums1, bial2_sub, colsums2)
R2  <- res[[1]]
M   <- res[[2]]

P   <- apply(M,1,function(x){return(fisher.test( matrix(x,ncol=2,byrow=TRUE) )$p.value ) 
       })


# Distance 
     if(length(colnames(bial1_sub))>1){
       names1 <- as.numeric(colnames(bial1_sub))
       names2 <- as.numeric(colnames(bial2_sub))
       d_dist <- apply(combinations,2,function(x){return( abs(names1[x[1]]-names2[x[2]]) )})
     }else{
       d_dist <- rep(NaN,dim(combinations)[2])
     }
      

res <- rbind(R2, P , d_dist)

 rownames(res)     <- c("R2","P","Distance")

# colnames(res)     <- apply(combinations,2,function(x){return(paste(x[1],"/",x[2],sep=""))})
# rownames(res)     <- c("d_raw","d_prime","r2","r","d_dist")
# colnames(res)     <- apply(combinations,2,function(x){return(paste(x[1],"/",x[2],sep=""))})

 reslist[[yy]]     <- res

}# End of iteration over pops

reslist           <- as.matrix(reslist)
rownames(reslist) <- paste("pop",1:npops)
colnames(reslist) <- "Linkdisequilibrium"

return(reslist)

}
