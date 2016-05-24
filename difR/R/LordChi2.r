# DIF: LORD'S CHI-SQUARED STATISTIC        

LordChi2<-function(mR,mF){
res<-NULL
mod<-as.character(ncol(mR))
model<-switch(mod,"2"="1PL","5"="2PL","6"="3PLc","9"="3PL")
for (item in 1:nrow(mR)){
if (model=="1PL"){
v1<-mR[item,1]
v2<-mF[item,1]
Sig1<-mR[item,2]^2
Sig2<-mF[item,2]^2
res[item]<-(v1-v2)^2/(Sig1+Sig2)
}
else{
if (model=="2PL" | model=="3PLc"){
v1<-mR[item,1:2]
v2<-mF[item,1:2]
Sig1<-cbind(c(mR[item,3]^2,mR[item,5]),c(mR[item,5],mR[item,4]^2))
Sig2<-cbind(c(mF[item,3]^2,mF[item,5]),c(mF[item,5],mF[item,4]^2))
}
else{
v1<-mR[item,1:3]
v2<-mF[item,1:3]
Sig1<-cbind(c(mR[item,4]^2,mR[item,7],mR[item,8]),c(mR[item,7],mR[item,5]^2,mR[item,9]),c(mR[item,8],mR[item,9],mR[item,6]^2))
Sig2<-cbind(c(mF[item,4]^2,mF[item,7],mF[item,8]),c(mF[item,7],mF[item,5]^2,mF[item,9]),c(mF[item,8],mF[item,9],mF[item,6]^2))
}
M<-solve(Sig1+Sig2)
V<-v1-v2
res[item]<-t(V)%*%M%*%V
}
}
return(as.numeric(res))
}




