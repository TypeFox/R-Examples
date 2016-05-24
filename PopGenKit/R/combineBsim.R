combineBsim <-
function(introduced, source, ndigit=3, subsample=0)
{

input=read.delim(file=introduced, sep='\r')
noext= gsub("[.].*$","",introduced)

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

titleline='File including one putative source and one simulated introduction'
header=input[1:(nloci),1]
header=as.matrix(header,ncol=1)
header=rbind(titleline,header)

popline='Pop'


#cette boucle pour diviser les données en populations
for (a in 1:npops)
{
filename=paste(noext,'_simul_rep', a, '.gen', sep='')

count=a
if (a==1) first=nloci+count+1
if (a>1) first=sum(popsizes[1:(count-1)])+nloci+count+1
last=sum(popsizes[1:count])+nloci+count
popdata=as.matrix(input[first:last,],ncol=1)

thispop=popsizes[count]

if (subsample>0) {
if (thispop>subsample) {
ech=sample(1:thispop, subsample, replace=F) 
popdata=as.matrix(popdata[ech,])
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

genotype=paste(allele1,allele2,sep='')

ind[x+1]=genotype
}
ind.all=rbind(ind.all,ind)
}
file.copy(from=source, to=filename)
write.table(popline, file=filename, append=T, quote=F, row.names=F, col.names=F, sep='\t')
write.table(ind.all, file=filename, append=T, quote=F, row.names=F, col.names=F, sep='\t')
}
}

