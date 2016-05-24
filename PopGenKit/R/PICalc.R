PICalc <-
function(datafile,ndigit=3)
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


PIC=rep(NA,times=nloci)


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

	PICterm1=0
	PICterm2=0
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
		
		PICterm1=(as.numeric(results[2*j,m])^2)+PICterm1
		}
	
		for (m in 1:(k-1))
		{

			for (n in (m+1):k)
			{
			PICterm2=(as.numeric(results[2*j,m])^2)*(as.numeric(results[2*j,n])^2)+PICterm2
			}
		}

	PIC[j]=1-PICterm1-(2*PICterm2)

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


write.table(table.results, file=filename, quote=F, row.names=F, col.names=T, sep='\t') 

locrow=c(1:nloci)

PICtable=cbind(locrow,PIC)


return('PIC'=PICtable)

}

