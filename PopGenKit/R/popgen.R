popgen <-
function(datafile,ndigit=3, freq.overall=T, freq.by.pop=T, genetic.stats=T, pairwise.fst=T)
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
#mettre les cellules NA vides pour lesthétisme !
for (r in 1:nrow(allpopresults))
{
for (c in 1:ncol(allpopresults))
{
if (is.na(allpopresults[r,c])==T) allpopresults[r,c]=''
}
} 

colnames(allpopresults)=col.name

filename2=paste(noext,'_Pops_freq.txt',sep='')
if (freq.by.pop==T) {
write.table(allpopresults, file=filename2, quote=F, row.names=F, col.names=T, sep='\t')
}

#ici toutes les variables à mettre dans le tableau (certaines requises pour alleles prives)

allpopn=NAallpopresults[,2]
allpopn=allpopn[!is.na(allpopn)]
allpopk=NAallpopresults[,4]
allpopk=allpopk[!is.na(allpopk)]
first.col=NAallpopresults[,1]
first.col=first.col[!is.na(first.col)]



#ici pour calculer Ho
if (genetic.stats==T) {
allpopHo=NULL
for (v in 1:npops)
{
popHo=rep(NA,times=nloci)

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
}}

allpopHo=cbind(allpopHo,popHo)
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


#ici pour calculer les alleles privés

allpopprivate=NULL
num.allpopk=as.numeric(allpopk)

for (v in 1:npops)
{
count=0
popno=v
otherpops=c(1:npops)
otherpops=otherpops[-v]
private.pop=rep(NA,times=nloci)
for (j in 1:nloci)
{ 
if (allpopk[((v-1)*nloci)+j]>0) {
vecalleleID=matpopresults[((2*j)-1),]
vecalleleID=vecalleleID[!is.na(vecalleleID)]
nbk=length(vecalleleID)

private=0
vec=rep(NA,times=length(otherpops))
for (i in 1:length(otherpops))
{
vec[i]=(2*(otherpops[i]-1)*nloci+(2*j))
}
freqvec=matpopresults[vec,]

for (k in 1:nbk)
{
alleleID=vecalleleID[k]
freqk=freqvec[,k]
freqk=as.numeric(freqk)
specific=as.numeric(matpopresults[(2*(v-1)*nloci+(2*j)),k])
if (specific>0) {
if (sum(freqk)==0) private=private+1 }
}
private.pop[j]=private
}
}
allpopprivate=cbind(allpopprivate,private.pop)
}


#ici pour calculer la richesse allélique
#avant de commencer les calculs, on doit savoir le minimum sample size pour chaque locus

minima=matrix(NA, nrow=npops, ncol=nloci)
for (v in 1:npops)
{
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
#if (sum(vec.missing)==popsizes[v]) minima[v,j]=NA else 
if (sum(vec.missing)==0) alleles2=alleles else alleles2=alleles[-vec.missing]
minima[v,j]=length(alleles2)/2 
}}
}

lowest=rep(NA,times=nloci)
for (i  in 1 : nloci)
{ lowest[i]=min(minima[,i], na.rm=T)  }

#calculer pour vrai

pop.richness=NULL

for (v in 1:npops)
{
popno=v
if (v==1) first=1
if (v>1) first=sum(popsizes2[1:(popno-1)])+1
last=sum(popsizes2[1:popno])
popdata=matinput[first:last,]
loc.richness=rep(NA, times=nloci)

for (j in 1:nloci)
{ 
if (allpopk[((v-1)*nloci)+j]>0) {
alleles=popdata[,j]
vec.missing=(which(alleles[]=='?'))
if (sum(vec.missing)==0) alleles2=alleles else alleles2=alleles[-vec.missing]
mat.alleles=matrix(alleles2, ncol=2, byrow=T)
sim=rep(NA,1000)
for (n in 1:1000)
{
 ech=sample(1:nrow(mat.alleles),lowest[j],replace=F)
resampled=mat.alleles[ech,]
sim[n]= length(unique(as.vector(resampled)))
}
loc.richness[j]=round(mean(sim), 2)
}}
pop.richness=cbind(pop.richness,loc.richness)

}




#créer le tableau de résultats
table.Ho=as.matrix(as.vector(allpopHo), ncol=1)
table.He=as.matrix(as.vector(allpopHe), ncol=1)
table.private=as.matrix(as.vector(allpopprivate), ncol=1)
table.richness=as.matrix(as.vector(pop.richness), ncol=1)



table.stats=cbind(allpopn,allpopk, table.richness,table.private, table.Ho,table.He)

for (p in 2:npops)
{ table.stats=rbind(table.stats[1:(((p-1)*nloci)+p-2),], NA, table.stats[(((p-1)*nloci)+p-1):nrow(table.stats),]) }

table.stats=rbind(NA,table.stats)
table.stats=cbind(first.col,table.stats)


col.name=rep('',times=ncol(table.stats))
col.name[1]= 'Locus#'
col.name[2]= 'n'
col.name[3]= 'k'
col.name[4]= 'Ar'
col.name[5]= 'priv'
col.name[6]= 'Ho'
col.name[7]= 'He'
colnames(table.stats)=col.name

filename3=paste(noext,'_Stats.txt',sep='')

#mettre les cellules NA vides pour lesthétisme !
for (r in 1:nrow(table.stats))
{
for (c in 1:ncol(table.stats))
{
if (is.na(table.stats[r,c])==T) table.stats[r,c]=''
}
} 

write.table(table.stats, file=filename3, quote=F, row.names=F, col.names=T, sep='\t')

}

#calculs de Fst

if (pairwise.fst==T) {

v1=c(1:(npops-1))
v2=c(2:npops)

pair=as.matrix(expand.grid(v1,v2,stringsAsFactors=T))

if (npops>2) {
remove=NULL
for (i in 1:nrow(pair))
{
un=pair[i,1]
deux=pair[i,2]
if (un>=deux) remove=c(remove,i)
}
pair2=pair[-remove,]
} else { pair2=pair   }

nbcomp=nrow(pair2)

all.loci.fst=NULL
all.sums=NULL
for (j in 1:nloci)
{
table.fst=matrix(NA,npops,npops)
table.sums=matrix(NA, npops, npops)
for (t in 1:nbcomp)
{
pop1=pair2[t,1]
pop2=pair2[t,2]

first1=(((2*pop1)-2)*nloci)+1
last1=2*pop1*nloci
pop1freqs=matpopresults[first1:last1,]
pop1counts=matpopcounts[first1:last1,]

first2=(((2*pop2)-2)*nloci)+1
last2=2*pop2*nloci
pop2freqs=matpopresults[first2:last2,]
pop2counts=matpopcounts[first2:last2,]

locus1freq=as.numeric(pop1freqs[j*2,])
locus2freq=as.numeric(pop2freqs[j*2,])

locus1counts=as.numeric(pop1counts[j*2,])
locus2counts=as.numeric(pop2counts[j*2,])

sum1=sum(locus1counts, na.rm=T)
sum2=sum(locus2counts, na.rm=T)

counts.all=locus1counts+locus2counts
sum.all=sum(counts.all, na.rm=T)
freq.all=counts.all/sum.all

He1=1-sum(locus1freq^2, na.rm=T)
He2=1-sum(locus2freq^2, na.rm=T)

Hs=((He1*sum1)+(He2*sum2))/sum.all
Ht=1-sum(freq.all^2, na.rm=T)
Fst=round(1-(Hs/Ht), 4)

table.fst[pop1,pop2]=Fst
table.sums[pop1,pop2]=sum.all

}

#ajouter les zeros sur la diagonale
for (i in 1:npops){
for (j in 1:npops) {
if (i==j) table.fst[i,j]=0 }}

all.loci.fst=rbind(all.loci.fst, table.fst)
all.sums=rbind(all.sums,table.sums)
}

#calculer multilocus Fst (moyenne pondérée)
multilocus.fst=matrix(NA, npops, npops)
for (t in 1:nbcomp)
{
pop1=pair2[t,1]
pop2=pair2[t,2]

weighted.freq=0
overall.sum=0
exclude.loci=0
for (j in 1:nloci) {
if (is.na(all.loci.fst[((npops*(j-1))+pop1),pop2])== T) 
{exclude.loci=exclude.loci+1}
else
{this.weighted.freq=all.loci.fst[((npops*(j-1))+pop1),pop2]*all.sums[((npops*(j-1))+pop1),pop2]
weighted.freq=weighted.freq+this.weighted.freq
overall.sum=overall.sum+all.sums[((npops*(j-1))+pop1),pop2] }}

multilocus.fst[pop1,pop2]=round(weighted.freq/overall.sum, 4)
if (exclude.loci>0) multilocus.fst[pop1,pop2]=paste(multilocus.fst[pop1,pop2], '!!', exclude.loci, sep='')

}

#ajouter les zeros sur la diagonale
for (i in 1:npops){
for (j in 1:npops) {
if (i==j) multilocus.fst[i,j]=0 }}



filename4=paste(noext,'_singlelocusFst.txt',sep='')

filename5=paste(noext,'_multilocusFst.txt',sep='')


for (t in 1:nbcomp)
{
pop1=pair2[t,1]
pop2=pair2[t,2]
for (j in 1:nloci) {
if (all.loci.fst[((npops*(j-1))+pop1),pop2]== 'NaN') all.loci.fst[((npops*(j-1))+pop1),pop2]= 'NoFst'
if (is.na(all.loci.fst[((npops*(j-1))+pop1),pop2])== T) all.loci.fst[((npops*(j-1))+pop1),pop2]= 'NoFst'
}}

#ajouter le nom des loci dans all.loci.fst
locus.name=NULL
for (j in 1:nloci) {
locus.name=c(locus.name,j,rep(NA, times=npops-1)) }

all.loci.fst=cbind(locus.name,all.loci.fst)

#mettre les cellules NA vides pour lesthétisme !
for (r in 1:nrow(all.loci.fst))
{for (c in 1:ncol(all.loci.fst)){
if (is.na(all.loci.fst[r,c])==T) all.loci.fst[r,c]='' }
} 

for (r in 1:nrow(multilocus.fst))
{for (c in 1:ncol(multilocus.fst)){
if (is.na(multilocus.fst[r,c])==T) multilocus.fst[r,c]='' }
} 

col.name=rep('',times=ncol(all.loci.fst))
col.name[1]= 'Locus#'
colnames(all.loci.fst)=col.name


write.table(all.loci.fst, file=filename4, quote=F, row.names=F, col.names=T, sep='\t')
write.table(multilocus.fst, file=filename5, quote=F, row.names=F, col.names=F, sep='\t')

}


return(list('matinput'=matinput,  'popsizes'=popsizes,  'matpopresults'=matpopresults, 'NAallpopresults'=NAallpopresults, 'matpopcounts'=matpopcounts, 'allpopk'=allpopk, 'all.loci.fst'=all.loci.fst, 'all.sums'=all.sums))
}

