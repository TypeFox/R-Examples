"addlines" <-
function(U,Y,col="green",lwd=1,lty=1) {
u1<-U[,1]; u2<-U[,2]
n<-dim(Y)[1]
for(i in 1:(n-1)){
for(j in (i+1):n){
if(!is.na(Y[i,j])) {
if(Y[i,j]!=0) {   segments(u1[i],u2[i],u1[j],u2[j],col=col,lwd=lwd) }
                     }
               }}
                              }

