SI3<-function(dd,pop,ind){
fls=as.factor(ind)
pop=as.factor(pop)
LL=levels(pop)
a=length(LL)
m=length(dd)
DT=data.frame(pop,fls,dd)
RRR=array(0,dim=c(1,a))
colnames(RRR)=LL
rownames(RRR)=c('SI3')
for (i in 1:a){
DD2=DT[DT$pop==LL[i],][,]
m2=length(DD2)
DDday=as.matrix(DD2[,3:m2])
DDpop=DD2[,1:2]
DDsub=data.frame(DDpop,DDday)
AA=as.matrix(DDsub$fls)
m3=length(levels(as.factor(AA)))
DDResult=array(0,dim=c(m3,m2))
colnames(DDResult)=colnames(DDsub)
for (j in 3:m2) {
BB=aggregate(DDsub[,j], list(DDsub$pop,DDsub$fls), sum)
DDResult[,1]=BB[,1]
DDResult[,2]=BB[,2]
DDResult[,j]=as.matrix(BB[,3])
}
print(LL[i])
DDsub2=t(DDResult[,3:m2])
DDsub3=unique(DDsub2)
colnames(DDsub3)=paste(LL[i],'-',levels(as.factor(AA)))
print(paste(LL[i],"Xip-","numbers of open flowers per day"))
print(DDsub3)
pairs(DDsub3)
DDsub4=cor(DDsub3)
print(paste(LL[i],"ri-"," all pairwise Pearson correlations coefficients (ri) of xit"))
print(DDsub4)
DDsub5=DDsub4[lower.tri(DDsub4)]
RRR[1,i]=mean(DDsub5)
}
print("SI3:the mean of ri")
(RRR)
list(Result=RRR)
}
