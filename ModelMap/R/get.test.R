
get.test<-function(	proportion.test,
			qdatafn=NULL,
			seed=NULL,
			folder=NULL,
			qdata.trainfn=paste(strsplit(qdatafn,split=".csv")[[1]],"_train.csv",sep=""),
			qdata.testfn=paste(strsplit(qdatafn,split=".csv")[[1]],"_test.csv",sep="")){
	
	
## Select dataset
if (is.null(qdatafn)){
	if(.Platform$OS.type=="windows"){
		## Adds to file filters to Cran R Filters table.
		Filters <- rbind(Filters,img=c("Imagine files (*.img)", "*.img"))
		Filters <- rbind(Filters,csv=c("Comma-delimited files (*.csv)", "*.csv"))

		qdatafn <- choose.files(caption="Select data file", filters = Filters["csv",], multi = FALSE)
		if(is.null(qdatafn)){stop("")}
	}else{stop("you must provide qdatafn")}
}

## check that qdata.trainfn and qdata.testfn are filenames

if(is.matrix(qdata.trainfn)==TRUE || is.data.frame(qdata.trainfn)==TRUE){
	stop("in the function get.test() 'qdata.trainfn' must be the filename for the new training dataset")}

if(is.matrix(qdata.trainfn)==TRUE || is.data.frame(qdata.trainfn)==TRUE){
	stop("in the function get.test() 'qdata.testfn' must be the filename for the new test dataset")}

## Check if file name is full path or basename
if(is.matrix(qdatafn)!=TRUE && is.data.frame(qdatafn)!=TRUE){
	if(identical(basename(qdatafn),qdatafn)){
		if(is.null(folder)){
			if(.Platform$OS.type=="windows"){
				folder<-choose.dir(default=getwd(), caption="Select directory")
			}else{
				folder<-getwd()}}
		qdatafn<-paste(folder,"/",qdatafn,sep="")
	}
}

## check if qdata.trainfn and qdata.testfn are full path or basename

if(identical(basename(qdata.trainfn),qdata.trainfn)){
	if(is.null(folder)){
		if(.Platform$OS.type=="windows"){
			folder<-choose.dir(default=getwd(), caption="Select directory")
		}else{
			folder<-getwd()}}
	qdata.trainfn<-paste(folder,"/",qdata.trainfn,sep="")
}

if(identical(basename(qdata.testfn),qdata.testfn)){
	if(is.null(folder)){
		if(.Platform$OS.type=="windows"){
			folder<-choose.dir(default=getwd(), caption="Select directory")
		}else{
			folder<-getwd()}}
	qdata.testfn<-paste(folder,"/",qdata.testfn,sep="")
}


## Read in data

if(is.matrix(qdatafn)==TRUE || is.data.frame(qdatafn)==TRUE){
	qdata<-qdatafn
}else{
	qdata<-read.table(file=qdatafn,sep=",",header=TRUE,check.names=FALSE)}

if(!is.null(seed)){
	set.seed(seed)}

train<-sample(1:nrow(qdata),round(nrow(qdata)*(1-proportion.test)))
qdata.train<-qdata[train,]
write.table(qdata.train,  file = qdata.trainfn, sep=",",append = FALSE,row.names=FALSE)
if(nrow(qdata.train)<nrow(qdata)){	
	qdata.test<-qdata[-train,]
	write.table(qdata.test,  file = qdata.testfn, sep=",",append = FALSE,row.names=FALSE)
}
}