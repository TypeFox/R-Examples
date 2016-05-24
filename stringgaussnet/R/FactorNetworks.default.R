FactorNetworks.default <-
function (x,Factor,method="SIMoNe",options=NULL, ...)
{
	optionValue<-getOption("stringsAsFactors")
	options(stringsAsFactors=F)
	Factor<-as.factor(Factor)
	Levels<-levels(Factor)
	if (length(Factor)!=nrow(x$DataExpression)){stop("Wrong length of factor.")}
	if (length(Levels)==1){stop("Please use a factor with more than one level.")}
	if(!(method %in% c("SIMoNe","WGCNA"))){stop("Wrong method used. Must be SIMoNe or WGCNA.")}
	FactorNets<-list()
	
	if (length(options)==0)
	{
		if (method=="SIMoNe"){options<-list(NEdges=NA,ClusterMethod="both",AddAnnotations=F,MartDataset="hsapiens_gene_ensembl")}
		if (method=="WGCNA"){options<-list(SoftThreshold=8,AThreshold=0.85,AddAnnotations=F,MartDataset="hsapiens_gene_ensembl")}
	}
	for (i in 1:length(options))
	{
		if(is.character(options[[i]])){options[[i]]<-paste("'",options[[i]],"'",sep="")}
	}
	options<-unlist(options)
	
	for (Level in Levels)
	{
		cat("Level: ",Level,"\n",sep="")
		LevelPos<-which(Factor==Level)
		FactorDEGeneExpr<-x
		FactorDEGeneExpr$DataExpression<-FactorDEGeneExpr$DataExpression[LevelPos,]
		FactorNets[[Level]]<-list()
		FactorNets[[Level]]$DEGeneExpr<-FactorDEGeneExpr
		FactorNets[[Level]]$Network<-LevelNet<-eval(parse(text=paste("get",method,"Net(FactorDEGeneExpr,",paste(paste(names(options),options,sep="="),collapse=","),")",sep="")))
	}
	class(FactorNets)<-"FactorNetworks"
	options(stringsAsFactors=optionValue)
	FactorNets
}
