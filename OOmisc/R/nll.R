nll <-
function(fit,y){
fit<-as.matrix(fit)
y<-as.matrix(y)
fity<-cbind(fit,y)
-(sum(log(1-fity[fity[,2]==0,1]))+sum(log(fity[fity[,2]==1,1])))
}
