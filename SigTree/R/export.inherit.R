

export.inherit <- function(tree, unsorted.pvalues, adjust=TRUE, side=1, 
	method="hommel", file="", test="Stouffers", frame=FALSE, branch="edge")
{
#This is what one can reference the graph to in order to see which of the internal nodes
#  contain which tips.  It exports a csv file that has the node number in column 1, its
#  pvalue in column 2, and then the names of the tips that are in its family in columns
#  3 and on.  The first (n.tips) rows will have only 1 member of their row - themselves
#  (because these are all tips).
#User function
#Argument: tree is a phylogenetic tree of class phylo.
#Argument: unsorted.pvalues (can be sorted as well) is a data frame with the tip
#  identifier in column 1 and the pvalue to be used in column 2
#Argument: adjust is a logical argument controlling p-value adjustment.  If FALSE, then
#  no p-value adjustment will take place.  If true, then side and method will determine
#  the flavor of the adjustment.
#Argument: side is either 1 or 2.  This corresponds to whether the p-values are 1-sided
#  or 2-sided.
#Argument: method is one of the methods found in p.adjust.methods and is the method to
#  be used in the p.adjust function. 
#Argument: file is the file path to export to.
#Argument: test is either "Stouffers" or "Fishers."  This is the p-value combination
#  to be used.
#Argument: frame is a logical of whether or not to return a data.frame object (in addition to creating the file)
#Argument: branch controls branch definition: "edge" or "node" are only options.  This
#   does not affect statistical methods, only how groups of tips are numbered in the returned object.

	#error checking
	if(class(tree)!="phylo")
	{
		return(cat("Error: Class of tree must be phylo.","\n"))
	}
    tree <- reorder(tree, order='cladewise')
	if(length(tree$tip.label)!=length(unsorted.pvalues[,1]))
	{
		return(cat("Error: There must be the same number of tip labels in tree as in unsorted.pvalues.","\n"))
	}else if(mean(sort(as.character(unique(unsorted.pvalues[,1])))==sort(as.character(unique(tree$tip.label))))!=1)
	{
		return(cat("Error: The tip labels in tree must have the same labels as the tip labels in unsorted.pvalues.","\n"))
	}
	if(min(unsorted.pvalues[,2])<0 | max(unsorted.pvalues[,2])>1)
	{
		return(cat("Error: P-values in unsorted.pvalues must be between greater than or equal to 0 and less than or equal to 1.","\n"))
	}else
	{
		f.one <- 1-100*.Machine$double.eps
		f.zero <- 100*.Machine$double.eps
		t <- unsorted.pvalues[,2]>=f.one
		unsorted.pvalues[t,2] <- f.one
		t <- unsorted.pvalues[,2]<=f.zero
		unsorted.pvalues[t,2] <- f.zero
	}
	if(adjust!=TRUE & adjust!=FALSE)
	{
		return(cat("Error: Value of adjust must be either TRUE or FALSE.","\n"))
	}
	if(side!=1 & side!=2)
	{
		return(cat("Error: Value of side must be either 1 or 2.","\n"))
	}
	if(is.numeric(method))
	{
		return(cat("Error: Value of method should be one of those found in p.adjust.methods"))
	}
	if(test!="Stouffers" & test!="Fishers")
	{
		return(cat("Error: Value of test must be either \"Stouffers\" or \"Fishers\".","\n"))
	}
	if(test=="Stouffers" & side==2)
	{
		cat("Caution: Stouffer's Method is designed for 1-sided p-values.", "\n\n")
	}
	if(test=="Fishers" & side==1)
	{
		cat("Caution: For Fisher's Method applied to one-tailed p-values, significance thresholds for small p-values (near 0) are more meaningful than for large p-values (near 1).", "\n\n")
	}
    if(frame!=TRUE & frame!=FALSE)
	{
		return(cat("Error: Value of frame must be either TRUE or FALSE.","\n"))
	}
	if(frame==FALSE & file=='')
	{
	    return(cat("Error: if frame=FALSE, file must be specified","\n"))
	}
	if(branch != "node" & branch != "edge")
	{
		return(cat("Error: Value of branch must be either 'node' or 'edge'.","\n"))
	}

	#Get n.tips,  n.total.nodes, index through functions.  Create inherit, a data frame
	#  with n.total.nodes rows and n.tips columns.  results is obtained through the
	#  result function.
	n.tips <- num.tips(tree)
	n.total.nodes <- num.total.nodes(tree)
	index <- index.matrix(tree)
	inherit <- data.frame(matrix(nrow=n.total.nodes, ncol=n.tips))
	results <- result(tree, unsorted.pvalues, test, adjust, side, method)
	#for loop that goes through each row in order to assign all of the tips belonging
	#  to each row (node).  tips is obtained by using the indices of all of the 1's from
	#  a given column of index to get all of the correct tip labels
	#  (from tree$tip.label).  length.tips is the length of tips.  tips is assigned to
	#  the correct row of inherit.  We use 1:length.tips to make sure that we assign the
	#  correct number of columns (and thus don't repeat).  For example, if the 2nd tips
	#  had 4 values, we would assign these 4 values to inherit[4,]'s first four columns
	#  and none after.   Each row of inherit represents a node and each column a tip.
	for(i in 1:n.total.nodes)
	{
		tips <- tree$tip.label[index[,i]==1]
		length.tips <- length(tips)
		inherit[i,1:length.tips] <- tips
	}
	#create a new data frame names new.inherit that has tree$tip.label, ie. the names of
	#  the tips in the first column and then the numbers from (n.tips+1) to n.total.nodes
	#  as the rest of the values of the first column.  We do this because we want the tip
	#  labels to be displayed instead of just 1:(n.tips+1) (what R knows them as).
	#  Afterwards, we want the numbers from (n.tips+1) to n.total.nodes.  This way, we can
	#  reference them to the node labels in the plot.  In the second column, we use
	#  results[,1], which is the pvalues corresponding to the first column.  In the third
	#  column and on, we have the inherit data frame we just created which contains all of
	#  the tips that each node has in its family.
	Branch <- as.character(c(tree$tip.label,(n.tips+1):n.total.nodes))
	new.inherit <- data.frame(Branch,results[,1], inherit)
	#label the first column of new.inherit as "Branch"; label the second column as
	#  "Stouffer's p-value" or "Fisher's p-value" depending on the value of test.
	if(test=="Stouffers")
	  {names(new.inherit)[2] <- "Stouffer's p-value"} else {
   	   names(new.inherit)[2] <- "Fisher's p-value"}
	
	# If branch='edge', need to re-number the branches
	if(branch=="edge")
	  {
     	Branch <- as.character(c(tree$tip.label,(n.tips+1):n.total.nodes))
        Branch[n.tips+1] <- 'root'
        f1 <- data.frame(Node=Branch[-c(1:(n.tips+1))])
		f2 <- data.frame(Node=tree$edge[,2], Edge=c(1:num.edges(tree)))
        f <- merge(f1,f2, sort=F)
		Branch[(n.tips+2):nrow(new.inherit)] <- f$Edge
	    new.inherit$Branch <- Branch
	   }
	
	
	
	#write new.inherit to the filepath.  na='' makes all of the NA's become blank spaces.
	if(!frame){write.csv(new.inherit, file, na='', row.names=FALSE)}
	# return data.frame object if requested
	if(frame){return(new.inherit)}
}


