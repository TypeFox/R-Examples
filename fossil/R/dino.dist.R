`dino.dist` <-
function(x,method=sorenson,type="dis") 
{
mat.names <- list(colnames(x),colnames(x))
sim<- matrix(,ncol(x),ncol(x),dimnames=mat.names)
for (i in 2:ncol(x)) {
for (j in 1:(i-1)) sim[i,j] <- method(x[,i],x[,j])
  }
if (type=="sim") typ <- "similarity" 
  else {typ <- "dissimilarity"
    sim<-1-sim}
sim.mat<-as.dist(sim)
attr(sim.mat, "method") <- method
attr(sim.mat, "type") <- typ
return(sim.mat)
}

