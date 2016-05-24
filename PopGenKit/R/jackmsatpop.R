jackmsatpop <-
function(datafile, ndigit=3, interval=1, nrepet=100, richness=F)
{

input=read.delim(file=datafile, sep='\r')
noext= gsub("[.].*$","",datafile)

popcount=rbind(which(input == 'pop', arr.ind=T),which(input=='Pop',arr.ind=T),which(input=='POP',arr.ind=T))
popcount=popcount[,1]
npops=NROW(popcount)

popcount=sort(popcount)
nloci=popcount[1]-1

popsizes=vector(length=npops)

for (h in 1:npops)
{
if (h<npops) popsizes[h]=popcount[h+1]-popcount[h]-1
if (h==npops) popsizes[h]=nrow(input)-popcount[h]
}

n.individuals=sum(popsizes)

datamatrix=matrix(NA, n.individuals, nloci*2)

input2=input[-popcount,]
input3=input2[-(1:nloci)]

#cette boucle pour chaque individu 
for (b in 1:n.individuals)
{
inddata=input3[b]
inddata=as.vector(inddata)
inddata=gsub( "[^[:alnum:]]", "", inddata)

#ceci pour savoir où commencent les allèles, au cas où nom de lind inclus alphanum
total=nchar(inddata)
diff=total-(ndigit*2*nloci)

#cette boucle pour chaque locus par individu
for (x in 1:nloci)
{
datamatrix[b,((2*x)-1)]=substr(inddata,diff+(2*ndigit*(x-1))+(ndigit-(ndigit-1)), diff+(2*ndigit*(x-1))+ndigit)
datamatrix[b,(2*x)]=substr(inddata,diff+(2*ndigit*(x-1))+ndigit+1, diff+(2*ndigit*(x-1))+(2*ndigit))
}

}

n.col=ncol(datamatrix)
npop=NROW(popsizes)
matrichness=NULL
matrichnessDev=NULL
for (v in 1:npop)
{
no=v
if (v==1) first=1
if (v>1) first=sum(popsizes[1:(no-1)])+1
last=sum(popsizes[1:no])
pop=datamatrix[first:last,]
compute=jackmsat(pop, interval, nrepet, richness)
#name=c('pop_',no, '.txt')
name=paste('pop_',no, '.txt', sep='')
Means2=t(compute$Means)
Stdev2=t(compute$Stdev)
if (richness==F) write.table(Means2, file=name, quote=F, sep='\t', row.names=F, col.names=F)
if (richness==F) write.table(Stdev2, file=name, append=T, quote=F, sep='\t', row.names=F, col.names=F)
if (richness==T) matrichness=cbind(matrichness, compute$Means[-1,-1])
if (richness==T) matrichnessDev=cbind(matrichnessDev, compute$Stdev[-1,-1])
}
popnames=c(1:npop)
if (richness==T) {
matrichness=rbind(popnames, matrichness)
matrichnessDev=rbind(popnames, matrichnessDev)
filename=paste('Richness for ', interval, 'ind.txt')
rows2=c('',compute$row.names)
matrichness2=t(matrichness)
matrichness2=rbind(rows2, matrichness2)
matrichnessDev2=t(matrichnessDev)
matrichnessDev2=rbind(rows2, matrichnessDev2)
write.table(matrichness2, file=filename, quote=F, sep='\t', row.names=F, col.names=F)
write.table(matrichnessDev2, file=filename, append=T, quote=F, sep='\t', row.names=F, col.names=F)
}}

