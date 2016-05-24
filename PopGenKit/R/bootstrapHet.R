bootstrapHet <-
function(datafile,ndigit=3, nrepet=1000)
{
input=read.table(file=datafile, sep='\t', colClasses='character')
noext= gsub("[.].*$","",datafile)

nloci=ncol(input)-2

#transformer pour trouver npops et popsizes
whereNbSamples=grep('NbSamples',input[,1])
npops=gsub( "[^0123456789]", "", input[whereNbSamples,1])
npops=as.numeric(npops)

whereSampleSize =grep('SampleSize',input[,1])
popsizes=as.numeric(c(gsub('[^0123456789]', "", input[whereSampleSize,1])))

#reconnaitre les noms des pops
whereSampleName =grep('SampleName',input[,1])
popnames=(c(gsub('SampleName=', "", input[whereSampleName,1])))

#creer une matrice nind x nloci seulement
matinput=input[,3:ncol(input)]
xsums=rep(NA,times=nrow(matinput))
for (i in 1:nrow(matinput)){ xsums[i]=sum(nchar(matinput[i,])) }
emptyrows=which(xsums==0)
matinput=matinput[-emptyrows,]

#determiner le nombre dalleles/locus pour formatter output
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


#creer le tableau de resultats
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

#trier les alleles en ordre croissant dans le output
for (j in 1:nloci)
{
ordre=order(results[(2*j)-1,])
results[(2*j)-1,]=results[(2*j)-1,ordre]
results[(2*j),]=results[(2*j),ordre]
}


#ajouter une colonne au debut avec le no de locus et le % de donnees manquantes
loc.col=NULL
missing.col=NULL
k.col=NULL
n.alleles.col=NULL

for (i in 1:nloci)
{
loc.col=c(loc.col,i,NA) 
missing.col=c(missing.col,missing[i],NA) 
k.col=c(k.col,nbk[i],NA)
n.alleles.col=c(n.alleles.col,n.alleles[i],NA)
}



#ici les analyses par population commencent

allpopresults=NULL
matpopresults=NULL
matpopcounts=NULL
popsizes2=2*popsizes

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

if (popresults[2*j,m]=='NaN') popresults[2*j,m]=0 
}

popnbk[j]= length(which(popresults[2*j,]>0))
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
}

NAallpopresults=allpopresults


#ici toutes les variables a mettre dans le tableau (certaines requises pour alleles prives)

allpopn=NAallpopresults[,2]
allpopn=allpopn[!is.na(allpopn)]
allpopk=NAallpopresults[,4]
allpopk=allpopk[!is.na(allpopk)]
first.col=NAallpopresults[,1]
first.col=first.col[!is.na(first.col)]



#ici pour calculer Ho
allpopHo=NULL
allpopSampN=NULL
for (v in 1:npops)
{
popHo=rep(NA,times=nloci)
popSampN=rep(NA, times=nloci)

popno=v
if (v==1) first=1
if (v>1) first=sum(popsizes2[1:(popno-1)])+1
last=sum(popsizes2[1:popno])
popdata=matinput[first:last,]

for (j in 1:nloci)
{ 
if (allpopk[((v-1)*nloci)+j]>0) {
alleles=popdata[,j]
vec.missing=(which(alleles[]=='?'))
if (sum(vec.missing)==0) alleles2=alleles else alleles2=alleles[-vec.missing]
true.n=length(alleles2)/2
count=0
for (n in 1:(true.n))
{
second.allele= alleles2[(2*n)]
if (alleles2[(2*n)-1]!=second.allele) count=count+1 
}

popHo[j]= round(count/true.n, 3)
popSampN[j]=true.n
}}

allpopHo=cbind(allpopHo,popHo)
allpopSampN=cbind(allpopSampN, popSampN)
}


#ici pour calculer He

allpopHe=NULL
for (v in 1:npops)
{
popHe=rep(NA,times=nloci)


popno=v
first=(((2*popno)-2)*nloci)+1
last=2*popno*nloci
freqdata=matpopresults[first:last,]

for (j in 1:nloci)
{ 
if (allpopk[((v-1)*nloci)+j]>0) {
popHe[j]=round(1-sum(as.numeric(freqdata[2*j,])^2, na.rm=T),3)
}}

allpopHe=cbind(allpopHe,popHe)
}


#cràö©er le tableau de ràö©sultats
table.Ho=as.matrix(as.vector(allpopHo), ncol=1)
table.He=as.matrix(as.vector(allpopHe), ncol=1)

MeanHo=rep(NA,times=npops)
MeanHe=rep(NA,times=npops)
for (v in 1:npops)
{
MeanHo[v]=round(sum(allpopHo[,v]*(allpopSampN[,v]/sum(allpopSampN[,v]))),3)
MeanHe[v]=round(sum(allpopHe[,v]*(allpopSampN[,v]/sum(allpopSampN[,v]))),3)
}

table.stats=cbind(table.Ho,table.He)

if (npops>1) {
for (p in 2:npops)
{ table.stats=rbind(table.stats[1:(((p-1)*nloci)+p-2),], NA, table.stats[(((p-1)*nloci)+p-1):nrow(table.stats),]) } }

table.stats=rbind(NA,table.stats)
table.stats=cbind(first.col,table.stats)


col.name=rep('',times=ncol(table.stats))
col.name[1]= 'Locus#'
col.name[2]= 'Ho'
col.name[3]= 'He'
colnames(table.stats)=col.name

filename3=paste(noext,'_bootHe.txt',sep='')

#mettre les cellules NA vides pour l¨v\u2260esthàö©tisme !
for (r in 1:nrow(table.stats))
{
for (c in 1:ncol(table.stats))
{
if (is.na(table.stats[r,c])==T) table.stats[r,c]=''
}
} 

write.table(table.stats, file=filename3, quote=F, row.names=F, col.names=T, sep='\t')


#return(list('matinput'=matinput,  'popsizes'=popsizes,  'matpopresults'=matpopresults, 'NAallpopresults'=NAallpopresults, 'matpopcounts'=matpopcounts, 'allpopk'=allpopk, 'all.loci.fst'=all.loci.fst, 'all.sums'=all.sums))


#ici les calculs de bootstrap commencent

#ici les analyses par population commencent


simHo=matrix(NA, nrow=npops, ncol=nrepet)
simHe=matrix(NA, nrow=npops, ncol=nrepet)

for (repn in 1:nrepet)
{
allpopresults=NULL
matpopresults=NULL
matpopcounts=NULL
popsizes2=2*popsizes

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

popdatacopy=popdata
ech=sample(1:(nrow(popdata)/2), (nrow(popdata)/2), replace=T)

for (x in 1:length(ech))
{
popdata[((2*x)-1),]=popdatacopy[((ech[x]*2)-1),]
popdata[(2*x),]=popdatacopy[(ech[x]*2),]
}

locsimHo=rep(NA,times=nloci)
for (j in 1:nloci)
{ 

if (allpopk[((v-1)*nloci)+j]>0) 
{
alleles=popdata[,j]
vec.missing=(which(alleles[]=='?'))
if (sum(vec.missing)==0) alleles2=alleles else alleles2=alleles[-vec.missing]
true.n=length(alleles2)/2
count=0
for (n in 1:(true.n))
{
second.allele= alleles2[(2*n)]
if (alleles2[(2*n)-1]!=second.allele) count=count+1 
}

locsimHo[j]= round(count/true.n, 3)
}


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

if (popresults[2*j,m]=='NaN') popresults[2*j,m]=0 
}

popnbk[j]= length(which(popresults[2*j,]>0))
}

matpopresults=rbind(matpopresults,popresults)

matpopcounts=rbind(matpopcounts,popcounts)


simHo[v,repn]=mean(locsimHo)
}


#ici pour calculer He


for (v in 1:npops)
{
locsimHe=rep(NA, times=nloci)
popno=v
first=(((2*popno)-2)*nloci)+1
last=2*popno*nloci
freqdata=matpopresults[first:last,]

for (j in 1:nloci)
{ 
if (allpopk[((v-1)*nloci)+j]>0) {
locsimHe[j]=round(1-sum(as.numeric(freqdata[2*j,])^2, na.rm=T),3)
}}
simHe[v,repn]=mean(locsimHe)
}
}

CIs=matrix(NA, nrow=npops, ncol=7)

for (j in 1:nrow(CIs))
{
CIs[j,1]=j
CIs[j,2]=MeanHo[j]
CIs[j,3]=round(quantile(simHo[j,], probs=0.025, type=8),3)
CIs[j,4]=round(quantile(simHo[j,], probs=0.975, type=8),3)
CIs[j,5]=MeanHe[j]
CIs[j,6]=round(quantile(simHe[j,], probs=0.025, type=8),3)
CIs[j,7]=round(quantile(simHe[j,], probs=0.975, type=8),3)
}

colnames(CIs)=c('Pop#', 'Mean_Ho', 'LBound_Ho', 'UBound_Ho', 'Mean_He', 'LBound_He', 'UBound_He')

write.table(CIs, file=filename3, quote=F, row.names=F, col.names=T, sep='\t', append=T)

#return(list('Hetotal'=simHe, 'Hototal'=simHo))
}

