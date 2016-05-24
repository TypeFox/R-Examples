summary.ShortPathSTRINGNet <-
function(object, ...)
{
	cat("Distances summary:\n")
	print(summary(object$Edges$Distance))
	NIntermediates <- unique(object$Edges$NIntermediates)
	DistSummaries <- matrix(NA,ncol=length(NIntermediates),nrow=5)
	colnames(DistSummaries) <- NIntermediates
	rownames(DistSummaries) <- c("Count","Min distance","Max distance","Mean distance","Median distance")
	for (NIntermediate in NIntermediates)
	{
		Edges <- object$Edges[which(object$Edges$NIntermediates==NIntermediate),]
		colname <- as.character(NIntermediate)
		DistSummaries["Count",colname] <- nrow(Edges)
		DistSummaries["Min distance",colname] <- min(Edges$Distance)
		DistSummaries["Max distance",colname] <- max(Edges$Distance)
		DistSummaries["Mean distance",colname] <- mean(Edges$Distance)
		DistSummaries["Median distance",colname] <- median(Edges$Distance)
	}
	cat("\nDistance for each number of intermediates:\n")
	print(DistSummaries)
}
