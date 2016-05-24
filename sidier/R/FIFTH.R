FIFTH <-
function(inputFile=NA,align=NA,saveFile=T,outname=paste(inputFile,"IndelDistanceFifthState.txt",sep="_"),addExtremes=F)
{
#require(ape)
if(is.na(inputFile)==TRUE&is.na(align[1])==TRUE) print("Error: Please, define either alignment or input file")
if(is.na(inputFile)==FALSE&is.na(align[1])==FALSE) print("Error: Please, define either alignment or input file")
if(is.na(inputFile)==FALSE&is.na(align[1])==TRUE)
align<-read.dna(format="fasta",file=inputFile) 				#read alignment
mat_alin<-as.matrix(as.character(align))					#alignment as matrix


gapMAL<-c()
for (i in 1:ncol(mat_alin))
if(paste(mat_alin[,i],collapse="")==paste(rep("-",nrow(mat_alin)),collapse=""))
gapMAL<-c(gapMAL,i)
if (is.null(gapMAL)==F) 
mat_alin<-mat_alin[,-c(gapMAL)] #removes positions full of gaps

seqNames<-row.names(mat_alin)

DISTANCE<-matrix(0,ncol=nrow(mat_alin),nrow=nrow(mat_alin))
colnames(DISTANCE)<-seqNames
row.names(DISTANCE)<-seqNames



for (i in 1:(nrow(mat_alin)-1))
for (j in (i+1):nrow(mat_alin))
	{

	if(i==1&j==(i+1)) print("Starting computation.......")
	if(i==1&j==(i+1)) print("...............0% completed")
	if(i==round(nrow(mat_alin)/4)&j==(i+1)) print("..............25% completed")
	if(i==round(nrow(mat_alin)/2)&j==(i+1)) print("..............50% completed")
	if(i==(nrow(mat_alin)-1)&j==(i+1)) print(".............100% completed")
	if(i==(nrow(mat_alin)-1)&j==(i+1)) print("FINISHED!!.................")

	sec1<-mat_alin[i,]
	sec2<-mat_alin[j,]

	recod1<-sec1
	recod2<-sec2

	for(l in 1:length(sec1))
		{
		if(sec1[l]=="-")
		recod1[l]<-0
		if(sec1[l]!="-")
		recod1[l]<-1
		if(sec2[l]=="-")
		recod2[l]<-0
		if(sec2[l]!="-")
		recod2[l]<-1
		}

	RECOD<-rbind(recod1,recod2)
	
	compactSIN00<-RECOD

	if(ncol(RECOD)>1)
		{
		borrar6<-c()
		for(l in 1:ncol(RECOD))
		if(paste(RECOD[,l],collapse="")==paste(rep(0,2),collapse=""))
		borrar6<-c(borrar6,l)
		if(is.null(borrar6)==F)
		compactSIN00<-as.matrix(RECOD[,-borrar6])
		}
	
	compactSIN11<-compactSIN00
	if(ncol(compactSIN00)>1)
		{
		borrar8<-c()
		for(l in 1:(ncol(compactSIN00)))
		if(paste(compactSIN00[,l],collapse="")==paste(c(1,1),collapse=""))
		borrar8<-c(borrar8,(l))
		if(is.null(borrar8)==F)
		compactSIN11<-as.matrix(compactSIN00[,-borrar8])
		}
	
	DISTANCE[i,j]<-ncol(compactSIN11)
	DISTANCE[j,i]<-ncol(compactSIN11)
	}

#outputs
if(saveFile==TRUE)
write.table(DISTANCE,file=outname)
DISTANCE
}
