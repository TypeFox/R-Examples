MK <- function(matrix,populations,biallelics=FALSE,subModel=TRUE){


if(biallelics){matrix <- matrix[,biallelics]} # only look at biallelic positions

count <- length(populations)
comb  <- combn(count,2)


out <- matrix(NA,dim(comb)[2],8)
colnames(out) <- c("Ps","Pn","Ds","Dn","NI","alpha","G","P")
nn <- paste("pop",comb[1,1],"/pop",comb[2,1],sep="")
if(dim(comb)[2]>1){ # more than 2 Populations
 for(xx in 2:dim(comb)[2]){
    m <- paste("pop",comb[1,xx],"/pop",comb[2,xx],sep="")
    nn <- c(nn,m)
 } 
}#END if
row.names(out) <- nn

 for(xx in 1:dim(comb)[2]){
   pop1 <- populations[[comb[1,xx]]]
   pop2 <- populations[[comb[2,xx]]] 
   
  if(length(pop1)!=0 & length(pop2)!=0){ 
   tt   <- mktest(matrix,pop1,pop2)
   out[xx,1] <- as.numeric(tt$Ps)
   out[xx,2] <- as.numeric(tt$Pn)
   out[xx,3] <- as.numeric(tt$Ds)
   out[xx,4] <- as.numeric(tt$Dn)
   out[xx,5] <- as.numeric(tt$NI)
   out[xx,6] <- as.numeric(tt$alphax)
   out[xx,7] <- as.numeric(tt$G)
   out[xx,8] <- as.numeric(tt$P)
  } 
 }

return(list(MC=out))

}

