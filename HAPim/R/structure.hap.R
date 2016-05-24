`structure.hap` <-
function(nbre.marq,nbre.all.marq){
nbre.type=2^nbre.marq
nbre.hap=prod(nbre.all.marq)
Bin.type=matrix(NA,nrow=nbre.type,ncol=nbre.marq)
Hap=matrix(NA,nrow=nbre.hap,ncol=nbre.marq)
Index=list()
Som.freq=matrix(0,nrow=nbre.type,ncol=nbre.type)

for(ma in (1:nbre.marq) ){
Bin.type[,ma]=rep( c(rep(0,2^(nbre.marq-ma)),rep(1,2^(nbre.marq-ma))), 2^(ma-1))
}


for(ma in (1:nbre.marq)){
rep.left=1
if(ma>1) rep.left=prod(nbre.all.marq[1:(ma-1)])
rep.right=1
if(ma<nbre.marq) rep.right=prod(nbre.all.marq[(ma+1):nbre.marq])
temp=NULL
for(all in 1:nbre.all.marq[ma]){
temp=c(temp,rep(all,rep.right))}
Hap[,ma]=rep(temp,rep.left)
}


Index[[1]]=matrix(NA,nrow=1,ncol=nbre.hap)
Index[[1]][1,]=1:nbre.hap
Index[[nbre.type]]=matrix(NA,nrow=nbre.hap,ncol=1)
Index[[nbre.type]][,1]=1:nbre.hap


for(ik in (2:(nbre.type-1))){
Bin=Bin.type[ik,]
One=(1:nbre.marq)[Bin==1]
One.hap=prod(nbre.all.marq[One])
Zero=(1:nbre.marq)[Bin==0]
Zero.hap=prod(nbre.all.marq[Zero])

Index[[ik]]=matrix(NA,nrow=One.hap,ncol=Zero.hap)

hap.colle=rep(NA,nbre.hap)
for( ih in 1:nbre.hap){
hap.colle[ih]=paste(Hap[ih,One],collapse="")}
hap.uni=unique(hap.colle)

for( ih in 1:length(hap.uni)){
     Index[[ik]][ih,]=(1:nbre.hap)[hap.colle==hap.uni[ih]]}


}

Som.freq[1,1]=1
for(ik in (2:nbre.type)){
Bin=Bin.type[ik,]
     temp=rep(0,nbre.type)
        for(ii in (1:nbre.type)){
        for(ma in 1:nbre.marq){
        temp[ii]=temp[ii]+2^(nbre.marq-ma)*max(0,(Bin[ma]-Bin.type[ii,ma])) }}
     temp=unique(temp)
Som.freq[ik,ik]=1
Ordk=sum(Bin.type[ik,])
    for( ii in (2:length(temp))){
    ii1=temp[ii]+1
    Ordii=sum(Bin.type[ii1,])
    Som.freq[ik,ii1]=(-1)^(Ordk-Ordii)}
}
Som.freq=t(Som.freq)

res=list("Bin.type"=Bin.type,"Index"=Index,"Som.freq"=Som.freq,"Hap"=Hap)
res
}

