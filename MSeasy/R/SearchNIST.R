#R function for connection between R and NIST mass spectral library search tool on Windows
# require NIST MS (Mass Spectral) Search Program Version 2.0 or higher
SearchNIST<-function(mspfile=NULL){

st<-strsplit(date(), " ")[[1]]
stBis<-strsplit(st[4], ":")[[1]]
Hour<-paste(stBis[1], stBis[2], stBis[3], sep="-")
Date<-paste(st[1], st[2], st[3], Hour, sep="_")
Mypath<-paste("output_SearchNIST", "_", "result", Date, sep="")
dir.create(Mypath)

if (is.null(mspfile)==TRUE) {
	rm(mspfile)
}
#System Check 
if ((Sys.info()["sysname"])!="Windows"){
	cat("Sorry, this function only works on Windows. Your OS is not Windows.\n")
	
}
else
{
#Search for MSSEARCH location should be c:/nist/mssearch
if (file.exists(file.path("C:", "Windows","win.ini"))==TRUE){
	#nistcheck<-pmatch("[NISTMS]", readLines("win.ini"))
	nistpath<-readLines(file.path("C:", "Windows","win.ini"))[pmatch("Path32", readLines(file.path("C:", "Windows","win.ini")))]
	if(is.na(nistpath)!=TRUE){
		nistpath<-file.path(unlist(strsplit(nistpath, split="="))[2])
		
	}
	if(is.na(nistpath)){
		nistpath<-readLines(file.path("C:", "Windows","win.ini"))[pmatch("Path16", readLines(file.path("C:", "Windows","win.ini")))]
		print("Your NIST MS program version is <2.0 or impossible to detect NIST")
		if(is.na(nistpath)){
			require(tcltk)
			nistpath<-tclvalue(tkchooseDirectory(initialdir=getwd(), title="Please, select the MS SEARCH directory, should be c:/nist/mssearch/"))
			print(nistpath)
		}
	}
	
	
}
if (file.exists(file.path("C:", "Windows","win.ini"))==FALSE){
	require(tcltk)
	nistpath<-tclvalue(tkchooseDirectory(initialdir=getwd(), title="Please, select the MS SEARCH directory, should be c:/nist/mssearch/"))
}
if(file.exists(file.path(nistpath,"AUTOIMP.MSD"))==TRUE){
	firstlocatorpath<-file.path(nistpath,"AUTOIMP.MSD")
	secondlocatorpath<-file.path(readLines(file.path(firstlocatorpath)))
}
if(file.exists(file.path(nistpath,"AUTOIMP.MSD"))==FALSE){
	#create AUTIMP.MSD file
	#secondlocatorpath<-paste(nistpath,"FILESPEC.FIL", sep="")
	secondlocatorpath<-"FILESPEC.FIL"
	zz<-file(file.path(nistpath,"AUTOIMP.MSD"))
	cat(secondlocatorpath, file = zz, sep = "\n")
	close(zz)
	firstlocatorpath<-file.path(nistpath,"AUTOIMP.MSD")
}
if(file.exists(file.path(nistpath,"SRCREADY.TXT"))==TRUE){
	unlink(file.path(nistpath,"SRCREADY.TXT"))
}
#create FILSPEC.FIL

	if(exists("mspfile")==FALSE){
		require(tcltk)
		mspfile<-tk_choose.files(default=paste(getwd(),"/*.msp",sep=""),caption="Please, select the MSP file", multi=FALSE, filters=matrix(c("your msp file","*.msp"),ncol=2))
		#zz<-file(paste(nistpath,"FILESPEC.FIL", sep=""), "w")
		zz<-file("FILESPEC.FIL", "w")
		cat(paste(file.path(mspfile),"Overwrite",sep=" "), file = zz, sep = "\n")
		cat(paste(23,62789), file = zz, sep = "\n")
		close(zz)
	}
		if(exists("mspfile")==TRUE){
		
		#zz<-file(paste(nistpath,"FILESPEC.FIL", sep=""), "w")
		zz<-file("FILESPEC.FIL", "w")
		cat(paste(file.path(mspfile),"Overwrite",sep=" "), file = zz, sep = "\n")
		cat(paste(23,62789), file = zz, sep = "\n")
		close(zz)
	}
code<-paste(file.path(nistpath,"nistms$.exe"), " /Instrument /par=2", sep="")
#system("c:/nist05/mssearch/nistms$.exe /Instrument /par=2")
system(code)
stopflag=FALSE
print("Please wait...")
Sys.sleep(3)

while(stopflag==FALSE){

 if(file.exists(file.path(nistpath,"SRCREADY.TXT"))==TRUE){
	stopflag=TRUE
	print("Done")
	
 }

 

}
 file.copy(file.path(nistpath,"SRCRESLT.TXT"),file.path(Mypath,"ResultsFromNIST.txt"))
print(paste("A file called ResultsFromNIST.txt has been generated in ", Mypath, sep=" "))
 #a parser for that file is needed
 #clean unused files from directories needed
}#end else system check

}#end function