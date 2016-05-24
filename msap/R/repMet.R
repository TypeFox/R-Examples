repMet <-function(dataM,groups,nDec){
	cat("Report of methylation levels \n")
	dataM[dataM==11]<-"u"
	dataM[dataM==10]<-"h"
	dataM[dataM==1]<-"i"
	dataM[dataM==0]<-"f"
	dataM <- split(dataM,groups)
	
	ns<-names(dataM)
	res<-matrix(ncol=length(ns), nrow=4)
	colnames(res)<-ns
	rownames(res)<-c("HPA+/MSP+ (Unmethylated)", "HPA+/MSP- (Hemimethylated)","HPA-/MSP+ (Internal cytosine methylation)","HPA-/MSP- (Full methylation or absence of target)")
	for(x in seq_along(dataM)){
	res[1,x]<- as.numeric(table(as.matrix(dataM[[x]]))["u"])
	res[2,x]<- as.numeric(table(as.matrix(dataM[[x]]))["h"])
	res[3,x]<- as.numeric(table(as.matrix(dataM[[x]]))["i"])
	res[4,x]<- as.numeric(table(as.matrix(dataM[[x]]))["f"])
	suma<- sum(res[,x])
	res[1,x] <- res[1,x]/suma
	res[2,x] <- res[2,x]/suma
	res[3,x] <- res[3,x]/suma
	res[4,x] <- res[4,x]/suma
	}
	print(res, digits=nDec)
  return(dataM)
}