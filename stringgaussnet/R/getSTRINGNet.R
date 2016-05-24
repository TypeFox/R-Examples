getSTRINGNet <-
function (DEGeneExpr,Identifier=0,NAdditionalNodes=NA,Species=9606,ConvertAliases=T,AddAnnotations=F,MartDataset="hsapiens_gene_ensembl")
{
	optionValue<-getOption("stringsAsFactors")
	options(stringsAsFactors=F)
	if (class(DEGeneExpr)!="DEGeneExpr")
	{
		stop("First argument is of a wrong class (different of DEGeneExpr)")
	}
	if (Identifier==0)
	{
		Identifiers=rownames(DEGeneExpr$DEGenesResults)
	}
	if(Identifier!=0)
	{
		Identifiers=DEGeneExpr$DEGenesResults[,Identifier]
	}
	if (is.null(Identifiers))
	{
		stop("Problem with genes identifiers")
	}
	
	cat("Retrieving identifiers...\n")
	IdentifiersString <- paste(Identifiers,collapse="%0D")
	BeginURI="http://string-db.org/api/tsv-no-header/resolveList?identifiers="
	BeginURI2="http://string-db.org/api/psi-mi-tab/interactionsList?identifiers="
	AdditionalURI="&additional_network_nodes="
	SpeciesURI="&species="
	IdsURI <- paste(BeginURI,IdentifiersString,SpeciesURI,Species,sep="") 
	if(nchar(IdsURI)>=8198)
	{
		InitialLength=length(Identifiers)
		IdentifiersString <- substr(IdentifiersString,1,8198-nchar(BeginURI)-nchar(SpeciesURI)-nchar(Species)-nchar(AdditionalURI)-4) 
		Identifiers=strsplit(IdentifiersString,"%0D")[[1]]
		Identifiers=Identifiers[1:(length(Identifiers)-1)]
		IdentifiersString <- paste(Identifiers,collapse="%0D")
		IdsURI <- paste(BeginURI,IdentifiersString,SpeciesURI,Species,sep="")
		warning("Too many gene identifiers were used (",InitialLength,"). Only the ",length(Identifiers)," first genes have been used.")
	}
	STRINGIds<- read.table(IdsURI,sep="\t",header=F,quote="",stringsAsFactors=F)[,1]
	
	if(!is.null(NAdditionalNodes)){if (is.na(NAdditionalNodes)){NAdditionalNodes=2*length(unique(STRINGIds))}}
	RemoveAdded=F
	if (is.null(NAdditionalNodes)){NAdditionalNodes=0;RemoveAdded=T}
	if (NAdditionalNodes>=10^4)
	{
		NAdditionalNodes=9999
		warning("Too many additional nodes. Reduced to ",NAdditionalNodes)
	}
	
	cat("Extracting STRING network from API...\n")
	InteractionsURI <- paste(BeginURI2,IdentifiersString,SpeciesURI,Species,AdditionalURI,NAdditionalNodes,sep="")
	Interactions=read.table(InteractionsURI,sep="\t",header=F,stringsAsFactors=F)
	
	FormatInteractions=Interactions[,c(3,4,ncol(Interactions))] 
	names(FormatInteractions)=c("node1","node2","Interactions")
	CorrespondancesScores=c(neighborhood="nscore",fusion="fscore",cooccurence="pscore",homology="hscore",coexpression="ascore",experimental="escore",knowledge="dscore",textmining="tscore",combined_score="score") 
	for (NomColonne in names(CorrespondancesScores))
	{
		FormatInteractions[[NomColonne]]=rep(0,nrow(FormatInteractions))
		Correspondance=CorrespondancesScores[NomColonne]
		for (i in 1:nrow(FormatInteractions))
		{
			DifferentsScores=strsplit(FormatInteractions$Interactions[i],"\\|")[[1]] 
			PosScore=grep(paste("^",Correspondance,":",sep=""),DifferentsScores)
			if (length(PosScore)==1)
			{
				FormatInteractions[i,NomColonne]=as.numeric(strsplit(DifferentsScores[PosScore],":")[[1]][2]) 
			}
		}
	}
	FormatInteractions$Interactions=NULL
	
	if (ConvertAliases)
	{
		FormatInteractions <- convertAliases(FormatInteractions)
	}
	
	if(RemoveAdded)
	{
		AllNodes<-unique(c(FormatInteractions$node1,FormatInteractions$node2))
		AddedNodes<-AllNodes[which(!(AllNodes %in% rownames(DEGeneExpr$DEGenesResults)))]
		FormatInteractions <- FormatInteractions[which(!(FormatInteractions$node1 %in% AddedNodes | FormatInteractions$node2 %in% AddedNodes)),]
	}
	
	GenesAnnotations <- NULL
	if (AddAnnotations)
	{
		cat("Adding annotations...\n")
		if (requireNamespace("biomaRt",quietly=TRUE)) {ensembl <- biomaRt::useMart("ensembl", dataset = MartDataset)} else {stop("biomaRt package must be installed to use this function")}
		GenesAnnotations <- getGenesInformations(unique(c(FormatInteractions$node1,FormatInteractions$node2)),ensembl)
	}
	
	STRINGNetwork <- STRINGNet(FormatInteractions,DEGeneExpr,GenesAnnotations)
	if (Identifier==0)
	{
		NewIdentifiers=rownames(STRINGNetwork$DEGenes)
	}
	if(Identifier!=0)
	{
		NewIdentifiers=STRINGNetwork$DEGenes[,Identifier]
	}
	MissingGenes<-Identifiers[which(!(Identifiers %in% NewIdentifiers))]
	if (length(MissingGenes)>0)
	{
		cat("Those genes are not included in the network: ",paste(MissingGenes,collapse=", "),"\n",sep="")
	}
	
	options(stringsAsFactors=optionValue)
	return(STRINGNetwork)
}
