summary.STRINGNet <-
function(object, ...)
{
	options(stringsAsFactors=F)
	AllNodes <- unique(c(object$Edges$node1,object$Edges$node2))
	InitNodes <- rownames(object$DEGenes)
	AddedNodes <- AllNodes[which(!(AllNodes %in% InitNodes))]
	
	AllInteractions <- object$Edges[which(object$Edges$Interaction!="combined_score"),]
	InitInteractions <- AllInteractions[which(AllInteractions$node1 %in% InitNodes & AllInteractions$node2 %in% InitNodes),]
	AddedInteractions <- AllInteractions[which(!(AllInteractions$node1 %in% InitNodes & AllInteractions$node2 %in% InitNodes)),]
	
	Interactions <- list("All interactions"=AllInteractions,"Interactions between initial nodes"=InitInteractions)
	if(length(AddedNodes)>0){Interactions[["Interactions with added nodes"]]=AddedInteractions}
	
	Summaries=list()
	for (Name in names(Interactions))
	{
		SummaryTable <- as.data.frame(t(as.matrix(table(Interactions[[Name]]$Interaction))))
		AddedRows <- as.data.frame(matrix(NA,ncol=ncol(SummaryTable),nrow=4))
		names(AddedRows) <- names(SummaryTable)
		SummaryTable <- rbind(SummaryTable,AddedRows)
		rownames(SummaryTable) <- c("Count","Min score","Max score","Mean score","Median score")
		for (InteractionType in names(SummaryTable))
		{
			Scores=Interactions[[Name]]$Score[which(Interactions[[Name]]$Interaction==InteractionType)]
			SummaryTable["Min score",InteractionType] <- min(Scores)
			SummaryTable["Max score",InteractionType] <- max(Scores)
			SummaryTable["Mean score",InteractionType] <- mean(Scores)
			SummaryTable["Median score",InteractionType] <- median(Scores)
		}
		Summaries[[Name]] <- SummaryTable
	}
	
	for (Name in names(Summaries))
	{
		cat(Name,":\n",sep="")
		print(Summaries[[Name]])
		if(Name != names(Summaries)[length(Summaries)]){cat("\n")}
	}
	options(stringsAsFactors=T)
}
