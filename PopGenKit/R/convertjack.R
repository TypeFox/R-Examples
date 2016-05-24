convertjack <-
function(datafile, ndigit=3, rm.loci)
{

input=read.delim(file=datafile, sep='\r')
noext= gsub("[.].*$","",datafile)

popcount=rbind(which(input == 'pop', arr.ind=T),which(input=='Pop',arr.ind=T),which(input=='POP',arr.ind=T))
popcount=popcount[,1]
npops=NROW(popcount)

popcount=sort(popcount)
nloci=popcount[1]-1

jacknloci=nloci-(length(rm.loci))

popsizes=vector(length=npops)

for (h in 1:npops)
{
if (h<npops) popsizes[h]=popcount[h+1]-popcount[h]-1
if (h==npops) popsizes[h]=nrow(input)-popcount[h]
}


#cette boucle pour diviser les donn<U+00E9>es en populations
outall=NULL
for (a in 1:npops)
{
outpop=NULL
count=a
if (a==1) first=nloci+count+1
if (a>1) first=sum(popsizes[1:(count-1)])+nloci+count+1
last=sum(popsizes[1:count])+nloci+count
popdata=(input[first:last,])

thispop=popsizes[count]

#cette boucle pour chaque individu de la population a
for (b in 1:thispop)
{
ind=matrix(0,2,jacknloci+2)
ind[1,1]=b
ind[1,2]=1
ind[2,1]= ''
ind[2,2]= ''
inddata=popdata[b]
inddata=as.vector(inddata)
inddata=gsub( "[^[:alnum:]]", "", inddata)
#ceci pour savoir o<U+00F9> commencent les all<U+00E8>les, au cas o<U+00F9> nom de l<U+2019>ind inclus alphanum
total=nchar(inddata)
diff=total-(ndigit*2*nloci)

#cette boucle pour chaque locus par individu
realloc=0
for (x in 1:nloci)
{
if (x %in% rm.loci ==F) {
realloc=realloc+1
ind[1,realloc+2]=substr(inddata,diff+(2*ndigit*(x-1))+(ndigit-(ndigit-1)), diff+(2*ndigit*(x-1))+ndigit)
ind[2,realloc+2]=substr(inddata,diff+(2*ndigit*(x-1))+ndigit+1, diff+(2*ndigit*(x-1))+(2*ndigit))
} }
#fini tous les locus de cet individu, donc ajouter aux donn<U+00E9>es de la pop
outpop=rbind(outpop,ind)

}

#ajouter le header de la pop
options(useFancyQuotes = "U+022")
name=paste('Pop-',a,sep='')
name=dQuote(name)
popline1=paste('SampleName=',name,sep='')
popline2=paste('SampleSize=',popsizes[a],sep='')
popline3=paste('SampleData= {')

popheader=matrix('',4,2+jacknloci)
popheader[1,1]=popline1
popheader[2,1]=popline2
popheader[3,1]=popline3

popend=matrix('',2,2+jacknloci)
lastline=paste('}')
popend[1,1]=lastline

completepop=rbind(popheader,outpop,popend)

#fini toute la pop donc ajouter aux res
outall=rbind(outall,completepop)
}

#ajouter le header du fichier
options(useFancyQuotes = "U+022")
title='ConvertedFile'
title=dQuote(title)
Secondline=paste('Title=', title, sep='')
Firstline='[Profile]'
Thirdline=paste('NbSamples=',npops,sep='')
Fourthline='Datatype=MICROSAT'
Fifthline='GenotypicData=1'
Sixthline='GameticPhase=0'
Seventhline='LocusSeparator=TAB'
Eighthline='RecessiveData=0'
if (ndigit==2) missing='00'
if (ndigit==3) missing='000'
missSymbol='?'
missSymbol=sQuote(missSymbol)
Ninthline=paste('MissingData=',missSymbol,sep='')
Tenthline='[Data]'
Eleventhline='[[Samples]]'

fileheader=matrix('',15,2+jacknloci)
fileheader[1,1]=Firstline
fileheader[3,1]=Secondline
fileheader[4,1]=Thirdline
fileheader[5,1]=Fourthline
fileheader[6,1]=Fifthline
fileheader[7,1]=Sixthline
fileheader[8,1]=Seventhline
fileheader[9,1]=Eighthline
fileheader[10,1]=Ninthline
fileheader[12,1]=Tenthline
fileheader[14,1]=Eleventhline

outall=rbind(fileheader,outall)

#Arlequin reconnait seulement ? comme missing, alors corriger

for (x in 1:nrow(outall))
   {
for (y in 1:ncol(outall))
  {
if  (outall[x,y]==missing) outall[x,y]='?'
}
}

filename=paste(noext,'_jack.arp',sep='')
write.table(outall, file=filename, quote=F, row.names=F, col.names=F, sep='\t')
return(outall)
}

