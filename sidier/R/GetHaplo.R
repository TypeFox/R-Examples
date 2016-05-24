GetHaplo <-
function(inputFile=NA,align=NA,saveFile=TRUE,outname="Haplotypes.txt",format="fasta",seqsNames=NA,silent=FALSE)
{
#require(ape)
if(is.na(inputFile)==TRUE&is.na(align[1])==TRUE) print("Error: Please, define either alignment or input file")
if(is.na(inputFile)==FALSE&is.na(align[1])==FALSE) print("Error: Please, define either alignment or input file")
if(is.na(inputFile)==FALSE&is.na(align[1])==TRUE)
align<-read.dna(file=inputFile,format="fasta")

mat_alin<-as.matrix(as.character(align))
FH<-FindHaplo(align=align,saveFile=saveFile,outname=outname)
Huniques<-c()
U<-unique(FH[,2])
for(i in 1:length(U))
Huniques<-c(Huniques,which(FH[,2]==U[i])[[1]])
out<-align[Huniques,]

if(is.na(seqsNames[1])==F)
{
if(seqsNames[1]=="Inf.Hap")
dimnames(out)[[1]]<-U
if(seqsNames[1]!="Inf.Hap")
dimnames(out)[[1]]<-seqsNames
}

if(saveFile==T)
write.dna(out,file=outname,format=format)

if(silent==FALSE)
print(paste(length(U)," different haplotypes found",ifelse(saveFile==T, {paste(", and saved in the file: \"",outname,"\"",sep="")},{paste(", but not saved in any file",sep="")}),sep=""),quote=F)


out

}
