# RAJU'S AREA
RajuZ<-function(mR,mF,signed=FALSE){
res<-NULL
mod<-as.character(ncol(mR))
model<-switch(mod,"2"="1PL","5"="2PL","6"="3PL")

if (model=="1PL" | signed){
indic.b<-switch(model,"1PL"=1,"2PL"=2,"3PL"=2)
indic.se<-switch(model,"1PL"=2,"2PL"=4,"3PL"=4)
bR<-mR[,indic.b]
se.bR<-mR[,indic.se]
bF<-mF[,indic.b]
se.bF<-mF[,indic.se]
H<-bF-bR
sH<-sqrt(se.bF^2+se.bR^2)
if (model!="3PL") res<-cbind(H=H,sH=sH,Z=H/sH)
else{
c<-mR[,6]
res<-cbind(H=(1-c)*H,sH=(1-c)*sH,Z=H/sH)
}
}
else{
for (i in 1:nrow(mR)){
aF<-mF[i,1]
aR<-mR[i,1]
bF<-mF[i,2]
bR<-mR[i,2]
se.aF<-mF[i,3]
se.bF<-mF[i,4]
cov.F<-mF[i,5]
se.aR<-mR[i,3]
se.bR<-mR[i,4]
cov.R<-mR[i,5]
Y<-aF*aR*(bF-bR)/(aF-aR)
if (exp(Y)==Inf){
H<--(bF-bR)
BF<--1
BR<-1
AF<-AR<-0
}
else{
H<-2*(aF-aR)*log(1+exp(Y))/(aF*aR)-(bF-bR)
BF<-1-2*exp(Y)/(1+exp(Y))
BR<--BF
AF<-2*(Y*exp(Y)/(1+exp(Y))-log(1+exp(Y)))/aF^2
AR<--aF^2*AF/aR^2
}
sH<-BF^2*se.bF^2+BR^2*se.bR^2+AF^2*se.aF^2+AR^2*se.aR^2+2*BF*AF*cov.F+2*BR*AR*cov.R
res<-rbind(res,c(H=H,sH=sqrt(sH),Z=H/sqrt(sH)))
}
if (model=="3PL"){
c<-mF[,6]
res[,1]<-res[,1]*(1-c)
res[,2]<-res[,2]*(1-c)
}
}
return(list(res=res,signed=signed))}


