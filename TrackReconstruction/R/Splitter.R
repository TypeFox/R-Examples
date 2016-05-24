
Splitter <- function(TagFile,Begin,End,RmL,Hz)#,AnimalID,Subsection="")
{
	if('Time' %in% colnames(TagFile))
	{
		DateTime<-paste(as.character(TagFile$Date),as.character(TagFile$Time))
		TagFile<-cbind(DateTime,TagFile[,!names(TagFile) %in% c("Time","Date")])
	}
	
	Begin<-as.character(Begin)
	End<-as.character(End)
	TagFile$DateTime<-as.character(TagFile$DateTime)
	outlist=vector("list",length(Begin))#initiate list to hold trips
	for(i in 1:length(Begin))
	{
		#subset data by Begin and End times
		outlist[[i]]<-TagFile[(which(TagFile$DateTime==Begin[i])-RmL*Hz/2):(which(TagFile$DateTime==End[i])+RmL*Hz/2),]
	}
#	for(i in 1:length(Begin))
#	{
#		num<-i
#		num=ifelse(num<10 & length(Begin)>10,paste("0",num,sep=""),num)
#			num=ifelse(num<100 & length(Begin)>100,paste("0",num,sep=""),num)
#			num=ifelse(num<1000 & length(Begin)>1000,paste("0",num,sep=""),num)
#			num=ifelse(num<10000 & length(Begin)>10000,paste("0",num,sep=""),num)
#		Nombre<-paste(AnimalID,Subsection,num,".csv",sep="")#Create a name for the file
#		write.table(outlist[[i]],Nombre,sep="\t",row.names=FALSE,quote=FALSE)#This will write this file in the working directory
#	}
return(outlist=outlist)
}
