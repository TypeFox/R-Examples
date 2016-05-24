SI2_onepop<-function(dd,ind){
fls=as.factor(ind)
m=length(dd[1,])
a=length(levels(fls))
popdd=array(0,dim=c(m,a))
colnames(popdd)=levels(fls)
rownames(popdd)=colnames(dd)
for (i in 1:m){
popdd[i,]=tapply(dd[,i],fls,sum)
}
f=array(0,dim=c(1,a))
colnames(f)=colnames(popdd)
for (i in 1:a){
f[,i]=length(popdd[,i][popdd[,i]!=0])
}
f=t(f)
popa=1/(a-1)
popaij=array(0,dim=c(a,a))
colnames(popaij)=levels(fls)
rownames(popaij)=levels(fls)
for (i in 1:a){
for (j in 1:a){
popaij[i,j]=length(popdd[,i][((popdd[,j]>0)&(popdd[,i]>0))])
}}
poptime=diag(popaij)
poptime=as.matrix(poptime)
ej=colSums(popaij)-poptime
ej=as.matrix(ej)
efj=ej/f
n=a
si=efj/(a-1)
colnames(si)=c("si")
print(si)
}