export.figtree <- function(tree, unsorted.pvalues, adjust=TRUE, side=1, 
	method="hommel", p.cutoffs=ifelse(rep(side==1, ifelse(side==1, 6, 3)),
	c(.01, .05, .1, .9, .95, .99), c(.01, .05, .1)), file="", 
	pal=ifelse(rep(side==1, ifelse(side==1, 1, length(p.cutoffs)+1)),
	"RdBu", rev(brewer.pal(length(p.cutoffs)+1,"Reds"))),
	test = "Stouffers", edge.label=TRUE, ignore.edge.length=FALSE,
	branch="edge")
{
#Creates a file that can be loaded in figtree.  The branches are colored in figtree and
#  the edges are annotated with the p-value (if edge.label is TRUE)
#User function
#Argument: tree is a phylogenetic tree of class phylo.
#Argument: unsorted.pvalues (can be sorted as well) is a data frame with the tip
#  identifier in column 1 and the pvalue to be used in column 2
#Argument: adjust is a logical argument controlling p-value adjustment.  If FALSE, then
#  no p-value adjustment will take place.  If true, then side and method will determine
#  the type of adjustment.
#Argument: side is either 1 or 2.  This corresponds to whether the p-values are 1-sided
#  or 2-sided.
#Argument: method is one of the methods found in p.adjust.methods and is the method to
#  be used in the p.adjust function.  
#Arguments: p.cutoffs is a vector of numbers, x, where 0<x<1 (not including 0 and 1) and
#  sorted in ascending order.  These are the p value cutoffs.  They affect how the edges
#  are colored. If pal is a divergent palette (such as RdBu; probably use a divergent
#  palette if p-values are 1-sided),   the smallest and the largest p values are colored
#  the darkest.  For example, if p.cutoffs<-c(.01, .05, .10, .90, .95, .99), then the
#  ranges are [0-.01), (.01-.05], ..., (.99-1.]  When using the "RdBu" pallete, the
#  range 0-.01 is the darkest red shade and .99-1 is the darkest blue shade.  A
#  sequential palette would probably be used when p-values are 2-sided.  The default 
#  argument for p.cutoffs is c(.01, .05, .10, .90, .95, .99) if side is 1 and
#  c(.01, .05, .1) if side is 2.
#Argument: file is the file path to export to.
#Argument: pal is a pallete from the package RColorBrewer.  Find a list of palettes by
#  display.brewer.all() or ?brewer.pal or brewer.pal.info.  pal may also be a vector of
#  colors.  The length needs to be one longer than the length of p.cutoffs (because
#  this is how many ranges of colors there are).  Colors in this vector need to be in
#  hexadecimal format.  For example, "#B2182B".  	#pal=ifelse(side==1, "RdBu", "Reds")
#  assigns the default argument of pal to be "RdBu" if side is 1 and "Reds" if side is
#  2.  This is because it is probably better to use a divergent palette if using 1-sided
#  p-values and a sequential palette if using 2-sided p-values.  The default
#  argument for pal is "RdBu" when side is 1 and "Reds" (but in reverse) when side is 2.
#  Note that the sequential palettes in RColorBrewer range from light to dark.  This 
#  means that the light colors correspond to low p-values.  This is why the default
#  palette when side is 2 is a reversed version of "Reds", so that the darker reds
#  correspond to the lower p-values.
#Argument: test is either "Stouffers" or "Fishers."  This is the p-value combination to
#  be used.
#Argument: edge.label is a logical argument.  When TRUE, the edges in FigTree will have
#  annotations with P-values corresponding to the p-value for the node they connect
#  that is closer to the root.
#Argument: ignore.edge.length only has an effect if the original tree had edge lengths
#  defined.  If it did, then setting this argument to FALSE ignores these edge lengths.
#Argument: branch controls branch definition: "edge" or "node" are only options.  This
#   does not affect statistical methods, only the colors used in edge coloring.

	#Error checking
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
	if(max(p.cutoffs)>=1 | min(p.cutoffs)<=0 | is.unsorted(p.cutoffs))
	{
		return(cat("Error: The values of p.cutoffs must be greater than 0, less than 1, and in ascending order.","\n"))
	}
	if(length(p.cutoffs)==0)
	{
		return(cat("Error: There must be at least one value of p.cutoffs.","\n"))
	}
	if(length(pal)==1 & is.character(pal))
	{
		if(!is.element(pal, rownames(brewer.pal.info)))
		{
			return(cat("Error: pal must be either a palette from RColorBrewer(to see a list: rownames(brewer.pal.info)) or a vector of colors.","\n"))
		}
	}else if (length(pal)!=length(p.cutoffs)+1)
	{
		return(cat("Error: the numbers of colors in pal must be one more than the number of values in p.cutoffs.","\n"))
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
	if(edge.label!=TRUE & edge.label!=FALSE)
	{
		return(cat("Error: Value of edge.label must be either TRUE or FALSE.","\n"))
	}
	if(ignore.edge.length!=TRUE & ignore.edge.length!=FALSE)
	{
		return(cat("Error: Value of ignore.edge.length must be either TRUE or FALSE.","\n"))
	}
	if(branch != "node" & branch != "edge")
	{
		return(cat("Error: Value of branch must be either 'node' or 'edge'.","\n"))
	}
	#Use tree to create tree4d, a tree of class phylobase
	tree4d <- phylo4d(tree)
	#Use tree4d to create tree4dext, a tree of class phyext.
	new("phylo4d_ext", tree4d)->tree4dext
	#Get n.total.nodes, n.tips, results, and edgecolor from functions.
	n.total.nodes <- num.total.nodes(tree)
	n.tips <- num.tips(tree)
	results <- result(tree, unsorted.pvalues, test, adjust, side, method)
	edgecolor <- edge.colors(tree, unsorted.pvalues, p.cutoffs, pal, test, adjust, 
		side, method, branch)
	#Make a new column (column 2) in results that gives us the edge to the left of the
	#  current node (ie. the row number).  Edge to the left refers to the edge
	#  connecting to the immediate ancestor.  We use the which.edge function to find
	#  this immediate ancestor.  Then label this column "Edge to Left."

	# The following line added 07.14.14 to avoid problems with new match.edge function from ape 3.1-3
	results[,2] <- match(1:n.total.nodes, tree$edge[,2])
		
	names(results)[2]<-"Edge to Left"
	#create a data frame called nodecolor filled with NA's and having n.total.nodes rows.
	nodecolor<-data.frame(matrix(NA, nrow=n.total.nodes))
				#it doesn't look like i need to include t here.
				#t <- !is.na(results[,2])
				#nodecolor[t,1] <- edgecolor[results[t,2]]
	#Assign edgecolor[results[,2]] to the first column of nodecolor.  results[,2] is the
	#  edge to the left of the current node (ie. row number)
	nodecolor[,1] <- edgecolor[results[,2]]
	#Create p.cutoffs.new, which is p.cutoffs with a 1 appended to the end.  Create
	#  n.cutoffs, which is the length of p.cutoffs.new.  cols is a vector of color
	#  values from RColorBrewer based on n.cutoffs and pal.
	p.cutoffs.new <- c(p.cutoffs, 1)
	n.cutoffs <- length(p.cutoffs.new)
	
	if(length(pal)==1)
	{	
			cols <- brewer.pal(n.cutoffs, pal)
			if(side==1){ cols[ceiling(n.cutoffs/2)] <- brewer.pal(7, "Greys")[2] }
	}else cols <- pal
#This is only for the branch that is to the left(away from the tips) of the root.  This
#  has not yet been handled, and will be black in figtree unless we do this for loop.
#  I am assigning the same color to this branch as to the two edges that are connected
#  to the root node.  This part is similar to the for loop in the edgecolor function.
	for(k in n.cutoffs:1)
	{
		if(results[n.tips+1,1]<=p.cutoffs.new[k])   nodecolor[n.tips+1,1]<-cols[k]
	}
	#to get p-value for node to the right of current node.  only happens if they set
	#  edge.label=TRUE  i.e., they want a p-value annotation for each edge.
	#If edge.label==TRUE, we do this.  Create a data frame called node.pval full of
	#  NA's and with n.total.nodes rows.  Set results[tree$edge[results[,2],1],1] to be
	#  the first column.  results[,2] is the edge to the left of the current row (node).
	#  For a given row, tree$edge[results[,2],1] is the node in the first column (i.e.,
	#  parent node) of tree$edge that corresponds to the edge from results[,2], i.e., it
	#  is the parent node of the current node(row).  Then
	#  results[tree$edge[results[,2],1],1] is the p-value of the direct
	#  ancestor/parent/left node of the given direct descendant/child/right node (the
	#  current row number).  Assign this to node.pval's first column.  We do all of 
	#  this because of how figtree handles its labels.  Figtree's annotations are 
	#  attached to the edges, not the nodes.  We want to be able to highlight an edge
	#   and find the p-value of the node to the left, ie the parent node.  Because of
	#   the way figtree's annotations work, we need to assign this p-value to the child
	#   node of the edge.  That is why we need to assign node.pval this way.
	#!#
	# BUT -- beginning in package version 1.2, we allow the definition of branch to be
	#     either 'node' (only option in previous package versions) or 'edge'.
	#     The edge coloring accounts for this now -- but the p-value annotation can now
	#     be edge-specific (as in FigTree).
	if(edge.label)	
	{
		node.pval <- data.frame(matrix(NA, nrow=n.total.nodes))	
		node.pval[,1] <- results[tree$edge[results[,2],1],1]
		#fill in p-value for the branch coming out of the root (away from the tips).
		#  The p-value here is the same as the p-value for the two children branches
		#  from the root (when branch="node").
		if(branch=="node"){node.pval[n.tips+1,1] <- results[n.tips+1,1]}
		if(branch=="edge"){node.pval[,1] <- results[,1]}
		names(node.pval)[1] <- "pvalue of node to right"
	}
	#Use gsub and regular expressions to get rid of the "#" in front of the hexadecimal
	#  color.  Place in column 2.
	gsub("#","",nodecolor[,1]) ->nodecolor[,2]
	#Use strtoi to convert the hexadecimal color to a decimal color.  Place in third
	#  column of nodecolor.
	strtoi(nodecolor[,2],16) -> nodecolor[,3]
	#Label the columns of nodecolor
	names(nodecolor)[1]<-"Hexadecimal color"
	names(nodecolor)[2]<-"Hex. color w/o #"
	names(nodecolor)[3]<-"Decimal color"
	#Put information in tdata(tree4dext), which is extra information about our tree.  If
	#  edge.label==FALSE, we only put nodecolor[,3] (the decimal color information) in
	#  tdata(tree4dext).  But if edge.label==TRUE, we first create a new data frame
	#  called df which contains nodecolor[,3] in the first column and node.pval[,1] in
	#  the second column and name those columns "nodecolor" and "pvalue."  Then we assign
	#  df to tdata(tree4dext).
	if(edge.label)	
	{
		df <- data.frame(nodecolor[,3], node.pval[,1])
		names(df)[1:2] <- c("nodecolor", "pvalue")
		tdata(tree4dext) <- df 
	} else {tdata(tree4dext) <- nodecolor[,3]}
	#This code with op.warn is to suppress the warning messages.  After we have run the
	#  line to get temptext, we reset the warn option to what it originally was.  Similar
	#  with op.scipen.  Save the original parameter as op.scipen.  It has to do with when
	#  R prints things in scientific notation.  We need it to not use scientific notation
	#  in order for the regular expressions to work.
	op.warn <- options()$warn
	options(warn=-1)
	op.scipen <- options()$scipen
	options(scipen=1000)
	#Create temptext, which is the text of write.nexus.simmap of tree4dext (in a vector
	#  format).  If edge.label is TRUE, then we use version 1.5.  This is because
	#  tdata(tree4dext) contains two files and version 1.5 handles this.  If edge.label is
	#  FALSE, we only have onle column of tdata(tree4dext) and we can use the older
	#  version.  The regular expression code is different based on the value of edge.label.
	if(edge.label) { capture.output(write.nexus.simmap(tree4dext, vers=1.5), file = NULL, 
		append = FALSE) -> temptext } else {capture.output(write.nexus.simmap(tree4dext), file = NULL, 
		append = FALSE) -> temptext}
	#Set the options back to what they were before we changed them.
	options(warn=op.warn, scipen=op.scipen)
	#Create length.temptext, which is simply the length of temptext.
	length(temptext) -> length.temptext
	#If edge.label==TRUE, we use this regular expression statement, which evaluates
	#  temptext[length.temptext-1].  This is the only part of temptext that we need to use
	#  regular expressions on.  We make it be in a format that figtree can read and use
	#  our color information.   We also need to test to see if tree4dext has edge lengths
	#  (if not, they are all NA's and show up as 0's).  If they do, we need to include this
	#  in the exported file for Figtree.  If not, we don't include the 0.  If
	#  ignore.edge.length is true, we go to the else statement that doesn't include the
	#  edge lengths.	
	if(edge.label)
	{	
		{if(hasEdgeLength(tree4dext) & ignore.edge.length==FALSE)
		{
			temptext[length.temptext-1] <- 
				gsub(":\\[&nodecolor=\\{([0-9]+)\\},pvalue=\\{([0-9\\.]+)\\}\\]([0-9\\.]+)", 
				"\\[&!color=#\\1, P-value=\"\\2\\\"]:\\3", temptext[length.temptext-1])
		}
		else
		{
			temptext[length.temptext-1] <- 
				gsub(":\\[&nodecolor=\\{([0-9]+)\\},pvalue=\\{([0-9\\.]+)\\}\\]([0-9\\.]+)", 
				"\\[&!color=#\\1, P-value=\"\\2\\\"]", temptext[length.temptext-1])	
		}}	
	} 
	#If edge.label==FALSE, then we go to this else statement, which uses regular
	#  expressions to change temptext[length.temptext-1] to be in a format that figtree
	#  can read and use the color information.  We also need to test to see if tree4dext
	#  has edge lengths (if not, they are all NA's and show up as 0's).  If they do, we
	#  need to include this in the exported file for Figtree.  If not, we don't include
	#  the 0.   If ignore.edge.length is true, we go to the else statement that doesn't
	#  include the edge lengths.
	else
	{
		{if(hasEdgeLength(tree4dext) & ignore.edge.length==FALSE)			
		{	
			temptext[length.temptext-1] <- gsub(":\\{([0-9\\.]+),*([0-9\\.]*)\\}", 
			"[&!color=#\\1]:\\2", temptext[length.temptext-1])
		}	
		else 
		{
			temptext[length.temptext-1] <- gsub(":\\{([0-9\\.]+),*([0-9\\.]*)\\}", 
			"[&!color=#\\1]", temptext[length.temptext-1])
		}}
			
	}
	write(temptext, file=file)
}
