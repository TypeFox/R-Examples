jackmsat <-
function(table, interval=1, nrepet=100, richness=F)
{
library(stats)
mat1=as.matrix(table)
x=nrow(mat1)
nloci=(ncol(mat1)/2)
result=NULL
resmeans=NULL
resstdev=NULL
x2=x/interval
x3=floor(x2)
if (richness==T) x3=1
for (i in 1:nloci)
{
locus=mat1[,((2*i)-1):(2*i)]
indnull=NULL
ind=nrow(locus)
for (q in 1:ind)
{
if (as.numeric(locus[q,1])<1) indnull=c(indnull,q)
}
if (NROW(indnull)>0) locus=locus[-indnull,]
ind2=nrow(locus)
means=rep(0,x3)
 for (j in 1:x3)
{
sim=rep(0,nrepet)
for (k in 1:nrepet)
{
if (ind2>=j*interval) ech=sample(1:ind2,j*interval,replace=F)
if (ind2>=j*interval) relocus=locus[ech,]
if (ind2>=j*interval) locvec=as.vector(relocus)
if (ind2>=j*interval) alleles=unique(locvec)
if (ind2>=j*interval) sim[k]=NROW(alleles)
}
result=cbind(result,sim)
}
means=colMeans(result)
stdev=sd(result)
}
resmeans=rbind(resmeans,means)
final=matrix(resmeans,nloci,x3,byrow=T)
resstdev=rbind(resstdev,stdev)
finalstdev=matrix(resstdev,nloci,x3,byrow=T)
final=round(final,2)
finalstdev=round(finalstdev,2)
col.names= ''
for (m in 1:x3)
{
col.names=c(col.names,m*interval)
}
row.names=NULL
for (n in 1:nloci)
{
row.names=c(row.names,n)
}
final=cbind(row.names,final)
final=rbind(col.names,final)
finalstdev=cbind(row.names,finalstdev)
finalstdev=rbind(col.names,finalstdev)
return(list(Means=final, Stdev=finalstdev, row.names=row.names, col.names=col.names))
}

