SI<-function(dd,pop){
pop=as.factor(pop)
m=length(colnames(dd))
a=length(levels(pop))
popdd=array(0,dim=c(m,a))
colnames(popdd)=levels(pop)
rownames(popdd)=colnames(dd)
for (i in 1:m){
popdd[i,]=tapply(dd[,i],pop,sum)
}
popa=1/(a-1)
popaij=array(0,dim=c(a,a))
colnames(popaij)=levels(pop)
rownames(popaij)=levels(pop)
popbij=array(0,dim=c(m,a))
for (i in 1:a){
for (j in 1:a){
popaij[i,j]=length(popdd[,i][((popdd[,j]>0)&(popdd[,i]>0))])
}}
poptime=diag(popaij)
popbij=array(0,dim=c(a,a))
colnames(popbij)=levels(pop)
rownames(popbij)=levels(pop)
for (i in 1:a){
for (j in 1:a){
popbij[i,j]=min(poptime[i],poptime[j])
}}
si=array(0,dim=c(1,a))
colnames(si)=levels(pop)
rownames(si)=c("Si")
for (i in 1:a){
for (j in 1:a){
si[i]=sum(popaij[i,j]/popbij[i,j])*popa
}}
print(paste("Result:S_mean=",round(mean(si),2),"S_sd=",round(sd(si),2)))
list(m=m,n=a,f=levels(pop),aij=popaij,bij=popbij,time=poptime,si=si)
}