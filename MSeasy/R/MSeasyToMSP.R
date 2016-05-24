#R function converting Output_peak.txt or Output_cluster.txt obtained by the function MS.clust of MSeasy into a MSP file for NIST mass spectral 

#The format MSP looks like
#Name: 752225-5899
#DB#: 7598
#Num Peaks: 77
#42 0.9856; 57 395.2356; ....
MSeasyToMSP<-function(filename="", outfilename="", cluster="", autosearch=FALSE){

st<-strsplit(date(), " ")[[1]]
stBis<-strsplit(st[4], ":")[[1]]
Hour<-paste(stBis[1], stBis[2], stBis[3], sep="-")
Date<-paste(st[1], st[2], st[3], Hour, sep="_")
Mypath<-paste("output_MStoMSP", "_", "result", Date, sep="")
dir.create(Mypath)


if ((filename!="")==TRUE){
filename<-filename
print(filename)
}
if ((filename!="")==FALSE||missing(filename)){
require("tcltk")
filename<-tk_choose.files(default=paste(getwd(),"/*.txt",sep=""),caption="Please, select the text file", multi=FALSE, filters=matrix(c("your text file","*.txt"),ncol=2))
}
if ((outfilename!="")==FALSE||missing(outfilename)){
outfile<-"ForNIST"
}
if ((outfilename!="")==TRUE){
outfile<-outfilename
}

d<-read.table(filename, header=FALSE, sep=" ") 
# in case you just want the spectrum of one cluster in one MSP file
if ((cluster!="")==TRUE){
beg<-which(d[1,]=="cluster")
if (length(beg)==0){
		beg<-which(d[1,]=="cluster_number")
		
		}
names<-d[1,]
if(length(cluster)>1){
d_temp<-NULL
	for(cl in 1:length(cluster)){
		d_temp<-rbind(d_temp,subset(d,d[,beg]==as.numeric(cluster[cl])))
	}
d<-d_temp
}
if(length(cluster)==1){

d<-subset(d,d[,beg]==as.numeric(cluster))
}
d<-rbind(names,d)

}



zz <- file(file.path(Mypath,paste(outfile,".msp",sep="")), "w")  # open an output file connection
for (j in 2:dim(d)[1]){
	mz<-list()
	k=1
	
	
#gestion the starting point depends on Output_peak.txt or Output_cluster.txt

start<-which(d[1,]==8)+1

flagname=0
if (length(start)==0) {
	start<-which(d[1,]=="8_first_mz")+1
	flagname=1
}

	for (i in start:dim(d)[2]){
		dname<-paste(d[j,2],"-",d[j,3],"-",d[j,1], sep="")
		if (flagname==1){
			dname<-paste("cluster_",d[j,1], sep="")
		}
		dDB<-d[j,1]
		if ((d[j,i]>0)==TRUE){
			mz[k]<-paste(d[1,i]," ",round(d[j,i],4),sep="")
			k=k+1
		}
		dPeak=k-1
		
	} 
	cat(" ", paste("Name: ", dname), paste("DB#: ", dDB), paste("Num Peaks: ",dPeak),  file = zz, sep = "\n")
	if (dPeak>20) {
		
		while (length(mz) >=21 ){
			
			cat(paste(mz[1:21]),"", file = zz, sep="; ")
			mz<-mz[22:length(mz)]
			cat(paste("\n"), file = zz)		
		}
		cat(paste(mz),"", file = zz, sep="; ")
		cat(paste("\n"), file = zz)
	}
	else{
	cat(paste(mz),"", file = zz, sep="; ")
	cat(paste("\n"), file = zz)
	}
	
} 
close(zz)
print(paste("A data file",outfile,".MSP has been generated in the folder:", Mypath, cat("\n")))
if(autosearch==FALSE){
print("Now, open your NIST MS Search Tool and launched the function SearchNIST on the .MSP file generated")
}
if(autosearch==TRUE){
print("Launch NIST search")
SearchNIST(mspfile=paste(file.path(Mypath,outfile),".msp",sep=""))
}
}
