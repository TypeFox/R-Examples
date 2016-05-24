areEventConfigurationsIdentical <- function(ephy, index1, index2){
	
	nodeset <- c(ephy$eventData[[index1]]$node, ephy$eventData[[index2]]$node);
	diffs <- sum(table(nodeset) == 1);
	return(diffs == 0);	
}
