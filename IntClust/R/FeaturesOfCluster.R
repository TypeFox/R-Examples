FeaturesOfCluster<-function(LeadCpds,Data,Threshold=1,Plot=TRUE,plottype="new",location=NULL){
	SubsetData=as.matrix(Data[which(rownames(Data)%in%LeadCpds),])
	
	Common=SubsetData%*%t(SubsetData)
	
	if(Threshold>length(LeadCpds)){
		stop("threshold is larger than the number of LeadCpds. This number can maximally be the number of LeadCpds.")
	}
	
	
	Features=colnames(SubsetData[,which(apply(SubsetData,2,sum)>=Threshold)])
	
	if(is.null(Features)){
		Features=""
		Plot=FALSE
	}
	
	#SharedAmongLeadCpds=SubsetData[,Features]
	
	if(Plot==TRUE){
		BinFeaturesPlot(LeadCpds=LeadCpds,OrderLab=rownames(Data),
				Features=as.character(Features),Data=Data,ColorLab=NULL,nrclusters=NULL,cols=NULL,name=c("Shared Features Among Selected Compounds"),
				margins=c(8.5,2.0,0.5,9.5),plottype=plottype,location=location)	
	}
	
	Out=list("Number_of_Common_Features"=Common,"SharedFeatures"=Features)
	
	return(Out)
}
