mismatch <- function(matrix_pol,populations,THETA){

npops    <- length(populations)
popnames <- paste("pop",1:npops)
if(npops!=length(THETA)){warning("You have to choose a theta value for every popualtion ")}

RE           <- matrix(NaN,3,npops)
rownames(RE) <- c("MDSD","MDG1","MDG2")
colnames(RE) <- popnames

for(xx in 1:npops){

if(length(populations[[xx]])==0){next;}
thetaT    <- THETA[xx]
popmatrix <- matrix_pol[populations[[xx]],,drop=FALSE]
poplength <- length(populations[[xx]])
if(poplength>1){seqpairs  <- combn(poplength,2)}else{seqpairs<-rbind(1,1)}# only one sequence
ncomb     <- choose(poplength,2)

res <- apply(seqpairs,2,function(x){

     seq1   <- popmatrix[x[1],]
     seq2   <- popmatrix[x[2],]
     id     <- which(!is.na(seq1) & !is.na(seq2) & seq1!=seq2)
     freq1  <- length(id)
     moment <- freq1 - thetaT
     freq2  <- moment*moment
     freq3  <- moment*moment*moment
     freq4  <- moment*moment*moment*moment
     
 return(c(freq2,freq3,freq4))
     
})
rownames(res) <- c("freq2","freq3","freq4")
E  <- rowSums(res)
s2 <- E["freq2"]
g1 <- E["freq3"]
g2 <- E["freq4"]

  value <- s2/(ncomb-1)
  if(!is.na(value)){
    s <- sqrt(s2/(ncomb-1))
	  RE["MDSD",xx] <- s
	}else{s <- FALSE}
	
	if(s){
    RE["MDG1",xx] <- g1*ncomb/((ncomb-2)*((ncomb-1)*s*s*s))
    RE["MDG2",xx] <- g2* (ncomb-1) * ncomb /((ncomb-3)*(ncomb-2)*(ncomb-1)*s*s*s*s)-(3*(ncomb-1)*(ncomb-1))/((ncomb-3)*(ncomb-2))
  }
}# END for over all POPS

return(RE)

}
