characteristicfunction <-
function(n,game){
game<-as.matrix(game)
coal<-coalitions(n)
matriz<-cbind(coal[[1]],coal[[2]],rep(0,2^n))
for (i in 2:nrow(matriz)){
aux<-which(matriz[i,]==1)
index<-which(matriz[i,n+1]==game[,1])
matriz[i,n+2]=game[index,2]
}
colnames(matriz)<-c(1:n,"Coalition","Cost")
rownames(matriz)<-rep(" ",2^n)
return(matriz)}
