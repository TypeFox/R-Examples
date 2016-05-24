MultiNetworks.default <-
function(x,Methods=c("STRING","SIMoNe","WGCNA"),STRINGOptions=NULL,SIMoNeOptions=NULL,WGCNAOptions=NULL,SelectInteractionsSTRING=NULL,STRINGThreshold=0,FilterShortPathOptions=NULL,FilterSIMoNeOptions=NULL,Factors=NULL, ...)
{
	optionValue<-getOption("stringsAsFactors")
	options(stringsAsFactors=F)
	if (class(x)=="DEGeneExpr"){x<-MultiDEGeneExpr(x)}
	MultiNets<-list()
	if(length(STRINGOptions)==0){STRINGOptions<-list(Identifier=0,NAdditionalNodes=NA,Species=9606,ConvertAliases=T,AddAnnotations=F,MartDataset="hsapiens_gene_ensembl")}
	if(length(SIMoNeOptions)==0){SIMoNeOptions<-list(NEdges=NA,ClusterMethod="both",AddAnnotations=F,MartDataset="hsapiens_gene_ensembl")}
	if(length(WGCNAOptions)==0){WGCNAOptions<-list(SoftThreshold=8,AThreshold=0.85,AddAnnotations=F,MartDataset="hsapiens_gene_ensembl")}
	for (ListName in names(x))
	{
		cat(ListName,"...\n\n",sep="")
		MultiNets[[ListName]]<-list()
		DEGene<-x[[ListName]]
		MultiNets[[ListName]]$DEGeneExpr<-DEGene
		for (method in Methods)
		{
			cat(method,"...\n\n",sep="")
			options<-eval(parse(text=paste(method,"Options",sep="")))
			for (i in 1:length(options)){if(is.character(options[[i]])){options[[i]]<-paste("'",options[[i]],"'",sep="")}}
			MultiNets[[ListName]][[method]]<-eval(parse(text=paste("get",method,"Net(DEGene,",paste(paste(names(options),options,sep="="),collapse=","),")",sep="")))
			if(method=="STRING")
			{
				if(length(SelectInteractionsSTRING)>0)
				{
					MultiNets[[ListName]][[method]]<-selectInteractionTypes(MultiNets[[ListName]][[method]],SelectInteractionsSTRING,STRINGThreshold)
				}
				newmethod<-"ShortPathSTRING"
				cat(newmethod,"...\n\n",sep="")
				MultiNets[[ListName]][[newmethod]]<-getShortestPaths(MultiNets[[ListName]][[method]])
				if(length(FilterShortPathOptions)>0)
				{
					cat("Filtering...\n\n")
					options<-FilterShortPathOptions
					for (i in 1:length(options)){if(is.character(options[[i]])){options[[i]]<-paste("'",options[[i]],"'",sep="")}}
					MultiNets[[ListName]][[newmethod]]<-eval(parse(text=paste("FilterEdges(MultiNets[[ListName]][[newmethod]],",paste(paste(names(options),options,sep="="),collapse=","),")",sep="")))
				}
			}
			if (method %in% c("SIMoNe","WGCNA") & !is.null(Factors))
			{
				newmethod<-paste(method,"FactorNetworks",sep="")
				cat(newmethod,"...\n\n",sep="")
				Factor<-Factors[rownames(DEGene$DataExpression)]
				if(length(Factor)!=nrow(DEGene$DataExpression)){stop("Wrong use of factors.")}
				MultiNets[[ListName]][[newmethod]]<-FactorNetworks(DEGene,Factor,method,eval(parse(text=paste(method,"Options",sep=""))))
			}
			if (method=="SIMoNe")
			{
				if (length(FilterSIMoNeOptions)>0)
				{
					cat("Filtering on global...\n\n")
					options<-FilterSIMoNeOptions
					for (i in 1:length(options)){if(is.character(options[[i]])){options[[i]]<-paste("'",options[[i]],"'",sep="")}}
					MultiNets[[ListName]][[method]]<-eval(parse(text=paste("FilterEdges(MultiNets[[ListName]][[method]],",paste(paste(names(options),options,sep="="),collapse=","),")",sep="")))
					newmethod<-paste(method,"FactorNetworks",sep="")
					if (!is.null(MultiNets[[ListName]][[newmethod]]))
					{
						cat("Filtering on ",newmethod,"...\n\n")
						MultiNets[[ListName]][[newmethod]]<-eval(parse(text=paste("FilterEdges(MultiNets[[ListName]][[newmethod]],",paste(paste(names(options),options,sep="="),collapse=","),")",sep="")))
					}
				}
			}
		}
	}
	class(MultiNets)<-"MultiNetworks"
	options(stringsAsFactors=optionValue)
	MultiNets
}
