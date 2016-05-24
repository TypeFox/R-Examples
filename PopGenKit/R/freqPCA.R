freqPCA <-
function(datafile, ndigit=3, default.pop=T, vecpops)
{
freq.overall=F
freq.by.pop=F

input=read.table(file=datafile, sep='\t', colClasses='character')
noext= gsub("[.].*$","",datafile)

nloci=ncol(input)-2

#transformer pour trouver npops et popsizes
whereNbSamples=grep('NbSamples',input[,1])
npops=gsub( "[^0123456789]", "", input[whereNbSamples,1])
npops=as.numeric(npops)

whereSampleSize =grep('SampleSize',input[,1])
popsizes=as.numeric(c(gsub('[^0123456789]', "", input[whereSampleSize,1])))

#reconnaître les noms des pops
whereSampleName =grep('SampleName',input[,1])
popnames=(c(gsub('SampleName=', "", input[whereSampleName,1])))

#créer une matrice nind x nloci seulement
matinput=input[,3:ncol(input)]
xsums=rep(NA,times=nrow(matinput))
for (i in 1:nrow(matinput)){ xsums[i]=sum(nchar(matinput[i,])) }
emptyrows=which(xsums==0)
matinput=matinput[-emptyrows,]

#déterminer le nombre dallèles/locus pour formatter output
kvec=vector(mode='numeric',nloci)
for (i in 1:nloci)
{
alleles=matinput[,i]
vec=unique(alleles)
vec=paste(vec,collapse='')
vec=gsub( "[^[:alnum:]]", "", vec)
k=nchar(vec)/ndigit
kvec[i]=k
}

MAX=max(kvec)

#créer le tableau de résultats
results=matrix(NA,2*nloci,MAX)

missing=rep(NA, times=nloci)
nbk=rep(NA, times=nloci)
n.alleles=rep(NA, times=nloci)
for (j in 1:nloci)
{ 
alleles=matinput[,j]
totaln=length(alleles)
vec=unique(alleles)
vec=paste(vec,collapse='')
vec=gsub( "[^[:alnum:]]", "", vec)
k=nchar(vec)/ndigit

sampsize=paste(alleles,collapse='')
sampsize=gsub( "[^[:alnum:]]", "", sampsize)
sampsize=(nchar(sampsize)/ndigit)

missingABS=length(grep('[?]',alleles))
missing[j]=round((100*(missingABS/totaln)),2)

nbk[j]=k

n.alleles[j]=sampsize/2

for (m in 1:k)
{
alleleID=substr(vec,(m*ndigit)-(ndigit-1),m*ndigit)
results[(2*j)-1,m]=alleleID
count=0
for (z in 1:length(alleles))
{
if (alleles[z]==alleleID) count=count+1
}
results[2*j,m]=round(count/sampsize,3)
}
}

#trier les allèles en ordre croissant dans le output
for (j in 1:nloci)
{
ordre=order(results[(2*j)-1,])
results[(2*j)-1,]=results[(2*j)-1,ordre]
results[(2*j),]=results[(2*j),ordre]
}

#ajouter une colonne au début avec le no de locus et le % de données manquantes
loc.col=NULL
missing.col=NULL
k.col=NULL
n.alleles.col=NULL

for (i in 1:nloci) {
loc.col=c(loc.col,i,NA) 
missing.col=c(missing.col,missing[i],NA) 
k.col=c(k.col,nbk[i],NA)
n.alleles.col=c(n.alleles.col,n.alleles[i],NA)
}

table.results=cbind(loc.col,n.alleles.col,missing.col, k.col, results)


#mettre les cellules NA vides pour lesthétisme !
for (r in 1:nrow(table.results))
{
for (c in 1:ncol(table.results))
{
if (is.na(table.results[r,c])==T) table.results[r,c]=''
}
} 

col.name=rep('',times=ncol(table.results))
col.name[1]= 'Locus#'
col.name[2]= 'n'
col.name[3]= 'Miss.%'
col.name[4]= 'k'
col.name[5]= 'Allele frequencies'
colnames(table.results)=col.name

filename=paste(noext,'_Overall_freq.txt',sep='')

if (freq.overall==T)  {
write.table(table.results, file=filename, quote=F, row.names=F, col.names=T, sep='\t') }

#ici les analyses par population commencent

allpopresults=NULL
matpopresults=NULL
matpopcounts=NULL
popsizes2=2*popsizes

PCA.table=NULL

for (v in 1:npops)
{
popmissing=rep(NA,times=nloci)
popnbk=rep(NA,times=nloci)
pop.n.alleles=rep(NA,times=nloci)

popresults=matrix(NA,2*nloci,MAX)
popcounts=matrix(NA,2*nloci,MAX)

popno=v
if (v==1) first=1
if (v>1) first=sum(popsizes2[1:(popno-1)])+1
last=sum(popsizes2[1:popno])
popdata=matinput[first:last,]

vecPCA=NULL

for (j in 1:nloci)
{ 
alleles=popdata[,j]
vecalleles=results[(2*j)-1,]
k=nbk[j]


sampsize=paste(alleles,collapse='')
sampsize=gsub( "[^0123456789]", "", sampsize)
sampsize=(nchar(sampsize)/ndigit)

missingABS=length(grep('[?]',alleles))
popmissing[j]=round((100*(missingABS/popsizes2[v])),2)

pop.n.alleles[j]=sampsize/2

for (m in 1:k)
{
alleleID=vecalleles[m]
count=0
for (z in 1:length(alleles))
{
if (alleles[z]==alleleID) count=count+1
}
popresults[(2*j),m]= round(count/sampsize,3)
popresults[(2*j)-1,]=vecalleles

popcounts[(2*j),m]= count
popcounts[(2*j)-1,]=vecalleles

if (popresults[2*j,m]=='NaN') popresults[2*j,m]='Miss!' 
}

popnbk[j]= length(which(popresults[2*j,]>0))

vecPCA2=as.numeric(popresults[(2*j),(1:(k-1))])
vecPCA=c(vecPCA,vecPCA2)
}

pop.missing.col=NULL
pop.k.col=NULL
pop.n.alleles.col=NULL

for (i in 1:nloci) {
pop.missing.col=c(pop.missing.col,popmissing[i],NA) 
pop.k.col=c(pop.k.col,popnbk[i],NA)
pop.n.alleles.col=c(pop.n.alleles.col,pop.n.alleles[i],NA)
}

table.popresults=cbind(loc.col,pop.n.alleles.col,pop.missing.col, pop.k.col, popresults)

blank.row=rep(NA,times=ncol(table.popresults))
blank.row[1]=popnames[v]
allpopresults=rbind(allpopresults,blank.row,table.popresults)

matpopresults=rbind(matpopresults,popresults)

matpopcounts=rbind(matpopcounts,popcounts)

PCA.table=rbind(PCA.table,vecPCA)
}

NAallpopresults=allpopresults
#mettre les cellules NA vides pour lesthétisme !
for (r in 1:nrow(allpopresults))
{
for (c in 1:ncol(allpopresults))
{
if (is.na(allpopresults[r,c])==T) allpopresults[r,c]=''
}
} 

colnames(allpopresults)=col.name

if (default.pop==F) {rownames(PCA.table)=vecpops} else {rownames(PCA.table)=popnames}

filename2=paste(noext,'_Pops_freq.txt',sep='')
if (freq.by.pop==T) {
write.table(allpopresults, file=filename2, quote=F, row.names=F, col.names=T, sep='\t')
}

filename3=paste(noext,'_FreqPCA.txt',sep='')
write.table(PCA.table, file=filename3, quote=F, row.names=T, col.names=F, sep='\t')

#return(list('matinput'=matinput,  'popsizes'=popsizes,  'matpopresults'=matpopresults, 'PCA'=PCA.table))
return(PCA.table)
}

