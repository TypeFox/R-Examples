LRmoments=function(p=c(0.5,0.5),kappaP=c(0,1,0),kappaD=c(1,0,0),log10=FALSE){
if(sum(p)<0|sum(p)>1) stop("Allele frequencies don't sum to 1")
if(abs(sum(kappaP)-1)>1e-6) stop("kappa-s should sum to 1")
if( kappaP[2]^2<4*kappaP[1]*kappaP[3]) stop("Violates kappa1^2>=4*kappa0*kappa2")
q=q012(p=p)
probHD=qkappa(kappa=kappaD,q=q) #Distribution under HD
probHP=qkappa(kappa=kappaP,q=q) #Distribution under HP
LR=probHP/probHD
if(log10) LR=log(LR,base=10)
muP=sum(LR*probHP,na.rm=TRUE)
muP.2=sum(LR^2*probHP,na.rm=TRUE)
SDP=sqrt(muP.2-muP^2)
skewP=sum(((LR-muP)/SDP)^3*probHP,na.rm=TRUE)
Plist=list(muP=muP,muP.2=muP.2,SDP=SDP,skewP=skewP)
muD=sum(LR*probHD)
muD.2=sum(LR^2*probHD)
SDD=sqrt(muD.2-muD^2)
skewD=sum(((LR-muD)/SDD)^3*probHD)
Dlist=list(muD=muD,muD.2=muD.2,SDD=SDD,skewD=skewD)
orderLR=order(LR)
pr.P=probHP[orderLR]
pr.D=probHD[orderLR]
LRdist=data.frame(LR=LR[orderLR],probHP=pr.P,probHD=pr.D)
aa=split(LRdist,LRdist$LR)
LRtable=NULL
for (i in 1:length(aa)){
 foo=apply(aa[[i]][,-1],2,sum)
 LRtable=rbind(LRtable,foo)
}
rownames(LRtable)=NULL
x=unique(LRdist[,1])
if (length(x)==dim(LRtable)[1])
 LRtable=data.frame(LR=unique(LRdist[,1]),LRtable)
else
 LRtable=LRdist
if(log10) cat("log10(LR):","\n")
list(moments=unlist(list(Plist,Dlist)),LRtable=LRtable)
}
