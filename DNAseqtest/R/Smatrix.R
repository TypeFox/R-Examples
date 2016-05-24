Smatrix <-
function(s,pix){
Pi<-diag(pix)
Sd<-matrix(0,4,4)
Sd[1,2:4]<-s[1:3]
Sd[2,3:4]<-s[4:5]
Sd[3,4]<-s[6]
for(i in 1:4){
Sd[,i]<-Sd[i,]
Sd[i,i]<-(-sum((Sd[i,-i])%*%(diag(Pi)[-i])))/diag(Pi)[i]
}
Sd
}
