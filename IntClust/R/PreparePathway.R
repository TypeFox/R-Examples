PreparePathway<-function(Object,GeneExpr,topG,sign){
	FoundGenes=NULL
	FoundComps=NULL
	
	FoundGenes=FindElement("Genes",Object)
	
	if(is.null(FoundGenes)|(is.list(FoundGenes) & length(FoundGenes) == 0)){
		FoundComps=FindElement("Compounds",Object)
		if(is.null(FoundComps)|(is.list(FoundComps) & length(FoundComps) == 0)){
			stop("Specify either the p-values of the genes or a selection of compounds to test for DE genes.")
		}
		
		pvalsgenes=list()
		FoundGenes=list()
		CompsP=list()
		TopDEP=list()
		for(i in 1:length(FoundComps)){
			LeadCpds=FoundComps[[i]]$LeadCpds
			CompsP[[i]]=FoundComps[[i]]
			names(CompsP)[[i]]=paste("Compounds_",i,sep="")
			if(is.null(LeadCpds)){
				stop("In the Compounds element, specify an element LeadCpds")
			} 
			
			group <- factor(ifelse(colnames(GeneExpr) %in% LeadCpds, 1, 0))
			
			if(class(GeneExpr)[1]=="ExpressionSet"){
				GeneExpr$LeadCmpds<-group		
				if (!requireNamespace("a4Base", quietly = TRUE)) {
					stop("a4Base needed for this function to work. Please install it.",
							call. = FALSE)
				}
					
				DElead <- a4Base::limmaTwoLevels(GeneExpr,"LeadCmpds")
				
				allDE <- a4Core::topTable(DElead, n = length(DElead@MArrayLM$genes$SYMBOL),sort.by="p")
				
				if(is.null(allDE$ID)){
					allDE$Genes <- rownames(allDE)
				}
				else
				{
					allDE$Genes=allDE$ID
				}
				
				if(is.null(topG)){
					topG=length(which(allDE$adj.P.Val<=sign))
				}
				TopDE <- allDE[1:topG, ]
				
				Genes <- list(TopDE,allDE)	
				names(Genes)<-c("TopDE","AllDE") 
			}
			else{				
				label.factor = factor(group)
				design = stats::model.matrix(~label.factor)
				fit = limma::lmFit(GeneExpr,design=design)
				fit = limma::eBayes(fit)
				allDE = limma::topTable(fit,coef=2,adjust="fdr",n=dim(GeneExpr)[1], sort.by="p")
				
				if(is.null(allDE$ID)){
					allDE$ID <- rownames(allDE)
				}
				
				if(is.null(topG)){
					topG=length(which(allDE$adj.P.Val<=sign))
				}
				
				
				TopDE<-allDE[1:topG,]
				
				Genes <- list(TopDE,allDE)	
				names(Genes)<-c("TopDE","AllDE") 
				
			}
			FoundGenes[[i]]=Genes
			TopDEP[[i]]=FoundGenes[[i]]
			names(TopDEP)[i]=paste("genes_",i,sep="")
			names(FoundGenes)=paste("Genes_",i,sep="")	
			pvalsgenes[[i]]=Genes$AllDE$P.Value
			names(pvalsgenes[[i]])=Genes$AllDE$ID
			names(pvalsgenes)[i]=paste("pvals_",i,sep="")
		}
	}	
	
	else{
		pvalsgenes=list()
		TopDEP=list()
		for(i in 1:length(FoundGenes)){
			#names(FoundGenes)[i]=paste("Genes_",i,sep="")
			TopDEP[[i]]=FoundGenes[[i]]
			names(TopDEP)[i]=paste("genes_",i,sep="")
			pvalsgenes[[i]]=FoundGenes[[i]]$AllDE$P.Value
			names(pvalsgenes[[i]])=FoundGenes[[i]]$AllDE$ID
			names(pvalsgenes)[i]=paste("pvals_",i,sep="")
		}
		
		FoundComps=FindElement("Compounds",Object)
		CompsP=list()
		for(i in 1:length(FoundComps)){
			if(is.null(FoundComps[[i]]$LeadCpds)){
				CompsP[[i]]="No LeadCpds specified"
			}
			else{
				CompsP[[i]]=FoundComps[[i]]
				names(CompsP)[[i]]=paste("Compounds_",i,sep="")
			}
		}
	}
	
	
	return(list(pvalsgenes=pvalsgenes,Compounds=CompsP,Genes=TopDEP))
}
