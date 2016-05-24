selectInteractionTypes <-
function(Network,InteractionTypes="All",Threshold=0)
{
	optionValue<-getOption("stringsAsFactors")
	options(stringsAsFactors=F)
	if (InteractionTypes[1]=="All"){InteractionTypes=unique(Network$Edges$Interaction)}
	KeepPos <- which(Network$Edges$Interaction %in% InteractionTypes & Network$Edges$Interaction!="combined_score" & Network$Edges$Score>=Threshold)
	if (length(KeepPos)==0)
	{
		stop("Did not find any interactions of interest")
	}
	SubNet <- Network$Edges[KeepPos,]
	ScoresSTRING <- matrix(NA,ncol=2+length(unique(SubNet$Interaction)),nrow=0)
	colnames(ScoresSTRING) <- c("node1","node2",unique(SubNet$Interaction))
	for (i in 1:nrow(SubNet))
	{
		node1 <- SubNet$node1[i]
		node2 <- SubNet$node2[i]
		RowPos <- which(ScoresSTRING[,"node1"]==node1 & ScoresSTRING[,"node2"]==node2)
		if (length(RowPos)==0)
		{
			NewRow <- c(node1,node2,rep(0,length(unique(SubNet$Interaction))))
			names(NewRow) <- colnames(ScoresSTRING)
			ScoresSTRING <- rbind (ScoresSTRING,NewRow)
			RowPos <- which(ScoresSTRING[,"node1"]==node1 & ScoresSTRING[,"node2"]==node2)
		}
		ScoresSTRING[RowPos,SubNet$Interaction[i]] <- SubNet$Score[i]
	}
	hscores <- rep(0,nrow(ScoresSTRING))
	if ("homology" %in% colnames(ScoresSTRING))
	{
		hscores <- ScoresSTRING[,"homology"]
	}
	rownames(ScoresSTRING) <- 1:nrow(ScoresSTRING)
	Scores <- as.data.frame(ScoresSTRING[,!(colnames(ScoresSTRING) %in% c("node1","node2","homology"))])
	for (j in 1:ncol(Scores)) {Scores[,j] <- as.numeric(Scores[,j])}
	Scores <- as.matrix(Scores)
	CombinedScores <- computeCombinedScores(Scores,hscores)
	
	SubNet<-Network
	SubNet$Edges<-SubNet$Edges[KeepPos,]
	AddedRows<-cbind(ScoresSTRING[,c("node1","node2")],Interaction=rep("combined_score",nrow(ScoresSTRING)),Score=CombinedScores)
	SubNet$Edges<-rbind(SubNet$Edges,AddedRows)
	SubNet$Edges$Score<-as.numeric(SubNet$Edges$Score)
	KeepNodes<-unique(c(SubNet$Edges$node1,SubNet$Edges$node2))
	SubNet$DEGenes<-SubNet$DEGenes[rownames(SubNet$DEGenes) %in% KeepNodes,]
	if (!is.null(SubNet$Annotations)){SubNet$Annotations<-SubNet$Annotations[rownames(SubNet$Annotations) %in% KeepNodes,]}
	
	options(stringsAsFactors=optionValue)
	return(SubNet)
}
