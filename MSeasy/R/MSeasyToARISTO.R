#R function converting Output_peak.txt or Output_cluster.txt from the function MS.clust of MSeasy into a file compatible with ARISTO webtool
#http://www.ionspectra.org/aristo/

#The format for Aristo batch search looks like
#>Name: 752225-5899
#42 0.9856
#57 395.2356
#..
MSeasyToARISTO<-function(filename="", outfilename="", cluster=""){

st<-strsplit(date(), " ")[[1]]
stBis<-strsplit(st[4], ":")[[1]]
Hour<-paste(stBis[1], stBis[2], stBis[3], sep="-")
Date<-paste(st[1], st[2], st[3], Hour, sep="_")
Mypath<-paste("output_MStoARISTO", "_", "result", Date, sep="")
dir.create(Mypath)

if ((filename!="")==TRUE){
filename<-filename
print(filename)
}
if ((filename!="")==FALSE||missing(filename)){
require("tcltk")
filename<-tk_choose.files(default=paste(getwd(),"/*.txt",sep=""),caption="Please, select a txt file", multi=FALSE, filters=matrix(c("your txt file","*.txt"),ncol=2))
}
if ((outfilename!="")==FALSE||missing(outfilename)){
outfile<-"ForARISTO"
}
if ((outfilename!="")==TRUE){
outfile<-outfilename
}

d<-read.table(filename, header=FALSE, sep=" ") 
#case you want just the spectras of one cluster in a ARISTO file
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


zz <- file(file.path(Mypath,paste(outfile,".txt",sep="")), "w")  # open an output file connection
for (j in 2:dim(d)[1]){
	mz<-list()
	k=1
#gestion of starting point 

start<-which(d[1,]==8)+1
flagname=0
if (length(start)==0){
 start<-which(d[1,]=="8_first_mz")+1
 flagname=1
 }
	for (i in start:dim(d)[2]){
		dname<-paste(d[j,2],"-",d[j,3],"-",d[j,1], sep="")
		if (flagname==1){
			dname<-paste("cluster_",d[j,1], sep="")
		}
		#dDB<-d[j,1]
		if ((d[j,i]>0)==TRUE){
			mz[k]<-paste(d[1,i]," ",round(d[j,i],4),sep="")
			k=k+1
		}
		dPeak=k-1
		
	} 
	cat(paste(">Name: ", dname),  file = zz, sep = "\n")
	cat(paste(mz),"", file = zz, sep="\n")
	#cat(paste("\n"), file = zz)
	
	
} 
close(zz)
print(paste("A data file", outfile,".txt has been generated in the folder:", Mypath, cat("\n")))
shell.exec("http://www.ionspectra.org/aristo/batchmode/")
print("Now, open www.ionspectra.org/aristo/batchmode/ file for submission")
}
