genLordChi2<-function(irtParam,nrFocal){
nrItems<-nrow(irtParam)/(nrFocal+1)
mat<-vector("list",nrFocal+1)
for (i in 1:(nrFocal+1)){
seq<-((i-1)*nrItems+1):(i*nrItems)
mat[[i]]<-irtParam[seq,]
}
mod<-as.character(ncol(irtParam))
model<-switch(mod,"2"="1PL","5"="2PL","6"="3PLc","9"="3PL")
C<-contrastMatrix(nrFocal,model)
nPar<-switch(model,"1PL"=1,"2PL"=2,"3PLc"=2,"3PL"=3)
res<-NULL
for (item in 1:nrItems){
Sig<-matrix(0,nPar*(nrFocal+1),nPar*(nrFocal+1))
v<-NULL
if (model=="1PL"){
for (gr in 1:(nrFocal+1)){
v<-c(v,mat[[gr]][item,1])
Sig[gr,gr]<-mat[[gr]][item,2]^2
}
}
else{
if (model=="2PL" | model=="3PLc"){
for (gr in 1:(nrFocal+1)){
v<-c(v,mat[[gr]][item,1:2])
si<-cbind(c(mat[[gr]][item,3]^2,mat[[gr]][item,5]),c(mat[[gr]][item,5],mat[[gr]][item,4]^2))
seq2<-((gr-1)*2+1):(gr*2)
Sig[seq2,seq2]<-si
}
}
else{
for (gr in 1:(nrFocal+1)){
v<-c(v,mat[[gr]][item,1:3])
si<-cbind(c(mat[[gr]][item,4]^2,mat[[gr]][item,7],mat[[gr]][item,8]),c(mat[[gr]][item,7],mat[[gr]][item,5]^2,mat[[gr]][item,9]),c(mat[[gr]][item,8],mat[[gr]][item,9],mat[[gr]][item,6]^2))
seq2<-((gr-1)*3+1):(gr*3)
Sig[seq2,seq2]<-si
}
}
}
Qj<-t(C%*%v)%*%solve(C%*%Sig%*%t(C))%*%(C%*%v)
res[item]<-Qj
}
return(res)
}



