FindHaplo <-
function(inputFile=NA,align=NA,saveFile=TRUE,outname="FindHaplo.txt")
{
#require(ape)
if(is.na(inputFile)==TRUE&is.na(align[1])==TRUE) print("Error: Please, define either alignment or input file")
if(is.na(inputFile)==FALSE&is.na(align[1])==FALSE) print("Error: Please, define either alignment or input file")
if(is.na(inputFile)==FALSE&is.na(align[1])==TRUE)
align<-read.dna(file=inputFile,format="fasta")

###
mat_alin<-as.matrix(as.character(align))
#seqs<-matrix(nrow=nrow(mat_alin),ncol=1)
#for(i in 1:nrow(mat_alin))
#seqs[i,]<-paste(mat_alin[i,],collapse="")

# To differentiate bad alignments on identical sequences:

SEQS<-list()

for (iseq in 1:nrow(mat_alin))
{
	if(length(which(mat_alin[iseq,]=="-"))>0)
	{
	SEQS[[iseq]]<-mat_alin[iseq,-c(which(mat_alin[iseq,]=="-"))]
	names(SEQS)[iseq]<-row.names(mat_alin)[iseq]
	}

	if(length(which(mat_alin[iseq,]=="-"))==0)
	{
	SEQS[[iseq]]<-mat_alin[iseq,]
	names(SEQS)[iseq]<-row.names(mat_alin)[iseq]
	}
}

seqs.lengths<-c()
for (iseq in 1:length(SEQS))
seqs.lengths<-c(seqs.lengths,length(SEQS[[iseq]]))

max.length<-max(seqs.lengths)

alin.new<-matrix("-",nrow=nrow(mat_alin),ncol=max.length)
for (iseq in 1:nrow(alin.new))
alin.new[iseq,1:seqs.lengths[iseq]]<-SEQS[[iseq]]

seqs<-rbind(paste(alin.new[1,],collapse=""),paste(alin.new[2,],collapse=""))
if(nrow(alin.new)>2)
for(i in 3:nrow(alin.new))
seqs<-rbind(seqs,paste(alin.new[i,],collapse=""))

###

UH<-unique(seqs)
NH<-length(UH)
namesH<-matrix(nrow=NH,ncol=1)
for(i in 1:NH)
namesH[i,]<-(paste("H",paste(rep(0,floor(log(NH,10))-floor(log(i,10))),collapse=""),i,sep=""))

ifelse(i<(10^floor(log(NH,10))), namesH[i,]<-(paste("H",paste(rep(0,floor(log(NH,10))-floor(log(i,10))),collapse=""),i,sep="")),
namesH[i,]<-(paste("H",i,sep=""))
)

out<-matrix(nrow=nrow(seqs),ncol=2)
colnames(out)<-c("Sequence.Name","Haplotype.Name")
out[,1]<-labels(align)

for(i in 1:NH)
out[which(seqs==UH[i,]),2]<-namesH[i,]

if(saveFile==T)
{
write.table(out,row.names=F,quote=F,file=outname)
}
out
}
