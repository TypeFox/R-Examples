fcons <-
function(indices,k){
samp<-c()
dim.ind<-dim(indices)
for(i in 1:(dim.ind[2])){
 if (length(unique(indices[,i]))<=k)
 samp<-c(samp,i)
}
indices.fil<-indices[,-samp]
indices.fil
}

