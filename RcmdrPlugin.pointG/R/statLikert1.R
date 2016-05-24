statLikert1<- function (X, y,tri=0) {
p <- ncol(X)
X <- as.matrix(X)

facto <- factor(y)
k <- length(levels(facto))
        
XXX <- matrix(0, nrow = k, ncol = p)
        
for (j in 1:p) {
    XXX[, j] <- tapply(X[, j], facto, mean, na.rm = TRUE)
        }
        
moyennes<-round(t(XXX),1)
moyennes<-as.data.frame(moyennes)
rownames(moyennes)<-colnames(X)
colnames(moyennes)<-levels(facto)
RC<-numeric(p)
p.value<-numeric(p)
for(j in 1:p){
ttt<-summary(aov(X[,j]~facto))[[1]]
RC[j]<-round(ttt[1,2]/ttt[2,2],3)
p.value[j]<-ttt[1,5]
}
moyennes<-data.frame(moyennes,RC,p.value)
if(tri==0){
return(moyennes)}
else{
return(moyennes[order(p.value),])}
}


