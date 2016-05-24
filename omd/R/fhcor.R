fhcor <-
function(indices,k){
cor.matrix<-cor(indices)
dim.cor<-dim(cor.matrix)
samp<-c()
for(i in 1:(dim.cor[1]-1)){
 for(j in (i+1):dim.cor[2]){
  if(abs(cor.matrix[i,j])>=k){
  samp<-c(samp,j);
  break}
}}
indices.fil<-indices[,-samp]
}

