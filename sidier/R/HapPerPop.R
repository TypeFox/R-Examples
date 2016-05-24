HapPerPop <-
function(inputFile=NA,sep=" ",input=NA,header=F,NameIniPopulations=NA, NameEndPopulations=NA,saveFile=T,Wname=NA,Iname=NA)
{

if(length(input)>1){
if(is.na(inputFile)==TRUE&is.na(input[1,1])==TRUE) print("Error: Please, define either input or input file")
if(is.na(inputFile)==FALSE&is.na(input[1,1])==FALSE) print("Error: Please, define either input or input file")
if(is.na(inputFile)==FALSE&is.na(input[1,1])==TRUE)
input<-read.table(inputFile,sep=sep,header=header)}

if(length(input)==1){
if(is.na(inputFile)==TRUE&is.na(input[1])==TRUE) print("Error: Please, define either input or input file")
if(is.na(inputFile)==FALSE&is.na(input[1])==FALSE) print("Error: Please, define either input or input file")
if(is.na(inputFile)==FALSE&is.na(input[1])==TRUE)
input<-read.table(inputFile,sep=sep,header=header)}

uhaplo<-sort(unique(input[,2]))
if(is.na(NameIniPopulations)==TRUE&is.na(NameEndPopulations)==TRUE)
pops<-input[,1]
if(is.na(NameIniPopulations)==FALSE&is.na(NameEndPopulations)==FALSE)
{
input[,1]<-substr(input[,1],NameIniPopulations,NameEndPopulations)
pops<-unique(input[,1])
#pops<-unique(substr(input[,1],NameIniPopulations,NameEndPopulations))
}

upops<-unique(pops)

matrizPresencia<-matrix(0,ncol=length(uhaplo),nrow=length(upops))
colnames(matrizPresencia)<-uhaplo
rownames(matrizPresencia)<-upops

for (i in 1:length(upops))
{
a<-input[which(upops[i]==input[,1]),2]
for(j in 1:length(a))
matrizPresencia[i,which(a[j]==colnames(matrizPresencia))]<-(1+matrizPresencia[i,which(a[j]==colnames(matrizPresencia))])
}

Pesos<-matrizPresencia
Interaccion<-matrizPresencia
Interaccion[which(matrizPresencia>0)]<-1

out<-list(c())
out[[1]]<-Pesos
out[[2]]<-Interaccion
names(out)<-c("Weighted","Interaction")

if(saveFile==T)
	{
		if(is.na(Wname))
			Wname<-"Weighted.txt"
		write.table(Pesos,file=Wname)

		if(is.na(Iname))
			Iname<-"Interaction.txt"
		write.table(Interaccion,file=Iname)
	}

out
}
