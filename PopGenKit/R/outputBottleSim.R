outputBottleSim <-
function(datafile, subsample=0)
{
ndigit=2
input=read.delim(file=datafile, sep='\r')
noext= gsub("[.].*$","",datafile)

load('alleles')


popcount=rbind(which(input == 'pop', arr.ind=T),which(input=='Pop',arr.ind=T),which(input=='POP',arr.ind=T))
popcount=popcount[,1]
npops=NROW(popcount)

popcount=sort(popcount)
nloci=popcount[1]-1

popsizes=vector(length=npops)

ref.alleles=c(11:72)

for (h in 1:npops)
{
if (h<npops) popsizes[h]=popcount[h+1]-popcount[h]-1
if (h==npops) popsizes[h]=nrow(input)-popcount[h]
}

titleline='File converted to original alleles from BottleSim'
header=input[1:(nloci),1]
header=as.matrix(header,ncol=1)
header=rbind(titleline,header)

filename=paste(noext,'_OkAllelesBsim.gen',sep='')

popline='Pop'

write.table(header, file=filename, quote=F, row.names=F, col.names=F, sep='\t')

input=as.matrix(input,ncol=1)

#cette boucle pour diviser les données en populations
for (a in 1:npops)
{
write.table(popline, file=filename, append=T, quote=F, row.names=F, col.names=F, sep='\t')

count=a
if (a==1) first=nloci+count+1
if (a>1) first=sum(popsizes[1:(count-1)])+nloci+count+1
last=sum(popsizes[1:count])+nloci+count
popdata=as.matrix(input[first:last,])

thispop=popsizes[count]

if (subsample>0) {
if (thispop>subsample) {
ech=sample(1:thispop, subsample, replace=F) 
popdata=popdata[ech,]
thispop=subsample }}

ind.all=NULL

#cette boucle pour chaque individu de la population a
for (b in 1:thispop)
{
ind=rep(NA, times=nloci+1)

inddata=popdata[b]
inddata=as.vector(inddata)

indID= gsub("[ ]","",inddata)
indID= gsub("[\t]|[,].*$","",indID)

indID=paste(indID,',',sep='')

ind[1]=indID

inddata=gsub( "[^[:alnum:]]", "", inddata)
#ceci pour savoir où commencent les allèles, au cas où nom de lind inclus alphanum
total=nchar(inddata)
diff=total-(ndigit*2*nloci)

#cette boucle pour chaque locus par individu
for (x in 1:nloci)
{
allele1=substr(inddata,diff+(2*ndigit*(x-1))+(ndigit-(ndigit-1)), diff+(2*ndigit*(x-1))+ndigit)
allele2=substr(inddata,diff+(2*ndigit*(x-1))+ndigit+1, diff+(2*ndigit*(x-1))+(2*ndigit))

if (x==1) vec.alleles=final.vec[1:kvec[1]]
if (x>1) vec.alleles=final.vec[(sum(kvec[1:(x-1)])+1):(sum(kvec[1:x]))]

for (k in 1:length(vec.alleles))
{
if (allele1==ref.alleles[k]) allele1=vec.alleles[k]
if (allele2==ref.alleles[k]) allele2=vec.alleles[k]
}

genotype=paste(allele1,allele2,sep='')

ind[x+1]=genotype
}
ind.all=rbind(ind.all,ind)
}
write.table(ind.all, file=filename, append=T, quote=F, row.names=F, col.names=F, sep='\t')
}
}

