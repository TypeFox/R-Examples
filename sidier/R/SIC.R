SIC <-
function(inputFile=NA,align=NA,saveFile=T,outnameDist=paste(inputFile,"IndelDistanceSIC.txt",sep="_"),outnameCode=paste(inputFile,"SIC_coding.txt",sep="_"),addExtremes=F)
{
#require(ape)
if(is.na(inputFile)==TRUE&is.na(align[1])==TRUE) print("Error: Please, define either alignment or input file")
if(is.na(inputFile)==FALSE&is.na(align[1])==FALSE) print("Error: Please, define either alignment or input file")
if(is.na(inputFile)==FALSE&is.na(align[1])==TRUE)
align<-read.dna(format="fasta",file=inputFile) 				#read alignment
mat_alin<-as.matrix(as.character(align))					#alignment as matrix

if(addExtremes==TRUE)
mat_alin<-cbind(rep("X",nrow(mat_alin)),mat_alin,rep("X",nrow(mat_alin)))

gapMAL<-c()
for (i in 1:ncol(mat_alin))
if(paste(mat_alin[,i],collapse="")==paste(rep("-",nrow(mat_alin)),collapse=""))
gapMAL<-c(gapMAL,i)
if (is.null(gapMAL)==F) 
mat_alin<-mat_alin[,-c(gapMAL)] #removes positions full of gaps

seqNames<-row.names(mat_alin)

	recod<-mat_alin
	if(
	length(which(recod[,ncol(recod)]=="-"))!=0|
	length(which(recod[,1]=="-"))!=0
	) stop(cat("\nSome sequences have gaps in extremes. Set 'addExtremes=TRUE' to analyse this dataset."))

	for(l2 in 1:nrow(mat_alin))
	for(l in 1:ncol(mat_alin))
		{
		if(mat_alin[l2,l]=="-")
		recod[l2,l]<-0
		if(mat_alin[l2,l]!="-")
		recod[l2,l]<-1
		}
	gapIni<-matrix(nrow=nrow(mat_alin),ncol=ncol(mat_alin))
	gapEnd<-matrix(nrow=nrow(mat_alin),ncol=ncol(mat_alin))

	for (i2 in 1:(ncol(recod)-1))
	for (i3 in 1:nrow(recod))
	{
	if(paste(recod[i3,i2],recod[i3,(i2+1)])=="1 0")
	gapIni[i3,i2]<-i2
	if(paste(recod[i3,i2],recod[i3,(i2+1)])=="0 1")
	gapEnd[i3,i2]<-i2
	}

	GAPS<-list()
	for (i4 in 1:nrow(gapIni))
	GAPS[[i4]]<-rbind(which(is.na(gapIni[i4,])==F),which(is.na(gapEnd[i4,])==F))

GAPSjoin<-cbind(GAPS[[1]],GAPS[[2]])
for(i5 in 3:length(GAPS))
GAPSjoin<-cbind(GAPSjoin,GAPS[[i5]])
GAPSunique<-unique(GAPSjoin,MARGIN=2)
GAPSunique<-GAPSunique[,order(GAPSunique[1,],GAPSunique[2,])] #sort by size

OUT<-rbind(GAPSunique,matrix(ncol=ncol(GAPSunique),nrow=nrow(mat_alin)))
row.names(OUT)<-c("Gap First Position","Gap Last Position",seqNames)
colnames(OUT)<-paste("Gap",c(1:ncol(OUT)),sep="")



for (i6 in 1:ncol(OUT))
for(i7 in 1:length(GAPS))
{
	coincide<-0
	if(length(GAPS[[i7]])>0)
	for (i8 in 1:ncol(GAPS[[i7]]))
	{
		if(paste(OUT[c(1:2),i6],collapse="-")==paste(GAPS[[i7]][,i8],collapse="-"))
		coincide<-1
	}
	OUT[(i7+2),i6]<-coincide
}


OUT<-as.data.frame(OUT)

for(i8 in 1:ncol(GAPSunique))
for(i9 in 1:ncol(GAPSunique))
if(i8!=i9)
if(GAPSunique[1,i8]<=GAPSunique[1,i9]&GAPSunique[2,i8]>=GAPSunique[2,i9])
{
OUT[(which(OUT[c(3:nrow(OUT)),i8]==1)+2),i9]<-"-"
}

DISTANCE<-matrix(0,ncol=nrow(mat_alin),nrow=nrow(mat_alin))
colnames(DISTANCE)<-seqNames
row.names(DISTANCE)<-seqNames

for(i10 in 3:(nrow(OUT)-1))
for(i11 in (i10+1):(nrow(OUT)))
	{
	seq1<-OUT[i10,]
	seq2<-OUT[i11,]
	
	Missing<-unique(c(which(seq1=="-"),which(seq2=="-"))) #Removes missings
	if(length(Missing)>0 & ncol(seq1)-length(Missing)==1)
	{
	seq1<-as.matrix(seq1[,-Missing])
	seq2<-as.matrix(seq2[,-Missing])
	}
	if(length(Missing)>0 & ncol(seq1)-length(Missing)>1)
	{
	seq1<-seq1[,-Missing]
	seq2<-seq2[,-Missing]
	}

	dis<-0
		for (i12 in 1:ncol(seq1))
		{
		if(seq1[,i12]!=seq2[,i12])
		dis<-(dis+1)
		}
DISTANCE[(i10-2),(i11-2)]<-dis
DISTANCE[(i11-2),(i10-2)]<-dis
	}

#outputs
#print(OUT)
#print(DISTANCE)
if(saveFile==TRUE)
{
write.table(OUT,outnameCode)
write.table(DISTANCE,outnameDist)
}
list(OUT,DISTANCE)
}
