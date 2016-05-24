plotSigTree <- function(tree, unsorted.pvalues, adjust=TRUE, side=1, 
	method="hommel", p.cutoffs=ifelse(rep(side==1, ifelse(side==1, 6, 3)),
	c(.01, .05, .1, .9, .95, .99), c(.01, .05, .1)),
	pal=ifelse(rep(side==1, ifelse(side==1, 1, length(p.cutoffs)+1)),
	"RdBu", rev(brewer.pal(length(p.cutoffs)+1,"Reds"))),
	test="Stouffers", branch.label=FALSE, tip.color=TRUE, edge.color=TRUE,
	tip.label.size=1, branch.label.size=1,  type="fan",
	use.edge.length=TRUE, edge.width=1, branch="edge", 
	root.edge=ifelse(type=="fan",FALSE,TRUE),
	branch.label.frame="none")
{
#Create plotSigTree, a function that plots our tree with optional edge colors and
#  tip colors.
# NOTE: Prior to package version 1.2, this called ape's plot.phylo function.
#       But beginning in version 1.2, it instead calls SigTree's plotphylo2 function.
#       This change allows full branch coloring (including the 'perpendicular-to-the-root' half-edges 
#       that were impossible to color differently in plot.phylo) and root edge coloring.
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
#  the darkest.  For example, if p.cutoffs = c(.01, .05, .10, .90, .95, .99), then
#  the ranges are [0-.01), (.01-.05], ..., (.99-1.]  When using the "RdBu" pallete, the
#  range 0-.01 is the darkest red shade and .99-1 is the darkest blue shade.  A
#  sequential palette would probably be used when p-values are 2-sided.  The default 
#  argument for p.cutoffs is c(.01, .05, .10, .90, .95, .99) if side is 1 and
#  c(.01, .05, .1) if side is 2.
#Argument: pal is a pallete from the package RColorBrewer.  Find a list of palettes by
#  display.brewer.all() or ?brewer.pal or brewer.pal.info.  pal may also be a vector of
#  colors.  The length needs to be one longer than the length of p.cutoffs (because this
#  is how many ranges of colors there are).  Colors in this vector need to be in
#  hexadecimal format.  For example, "#B2182B".  	#pal=ifelse(side==1, "RdBu", "Reds")
#  assigns the default argument of pal to be "RdBu" if side is 1 and "Reds" if side is 2.
#  This is because it is probably better to use a divergent palette if using 1-sided
#  p-values and a sequential palette if using 2-sided p-values.   The default
#  argument for pal is "RdBu" when side is 1 and "Reds" (but in reverse) when side is 2.
#  Note that the sequential palettes in RColorBrewer range from light to dark.  This 
#  means that the light colors correspond to low p-values.  This is why the default
#  palette when side is 2 is a reversed version of "Reds", so that the darker reds
#  correspond to the lower p-values.
#Argument: test is either "Stouffers" or "Fishers."  This is the p-value combination to
#  be used.
#Argument: branch.label is a logical parameter that displays numerical labels for the
#  branches (node or edge, depending on the branch argument) when TRUE
#  and doesn't when FALSE.  These labels can be used to cross-reference to the file obtained
#  from the export.inherit function.
#Argument: tip.color is a logical argument that enables the tips to be colored when TRUE and
#  not be colored when FALSE.
#Argument: edge.color is a logical argument that enables the edges to be colored when TRUE
#  and not be colored when FALSE.
#Argument: tip.label.size is a numerical argument that controls the (cex) size of the text
#  of the tip labels.
#Argument: branch.label.size is a numerical argument that controls the (cex) size of the text
#  of the branch labels (edge or node, depending on the branch argument).
#Argument: type controls which type of plot will be produced.  Available options are
#  "phylogram," "cladogram," "fan," "unrooted," and "radial."
#Argument: use.edge.length is a logical argument.  If TRUE, then the plot uses the defined
#	edge lengths as usual.  But if FALSE, then the plot ignores these edge lengths.  This
#	can be useful to produce a more uniformly-spaced tree.
#Argument: edge.width controls width of plotted edges. This is passed to plotphylo2
#
#Argument: branch controls branch definition: "edge" or "node" are only options.  This
#   does not affect statistical methods, only the colors used in edge coloring.
#
#Argument: root.edge is a logical argument. If TRUE, then the root edge is plotted.
#
#Argument: branch.label.frame is a character argument specifying how to frame branch labels
#    (passed to nodelabels or edgelabels depending on branch definition) -- must be "none", "circ", or "rect"
	
	#Error checking.
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
	if(length(p.cutoffs)==0)
	{
		return(cat("Error: There must be at least one value of p.cutoffs.","\n"))
	}
	if(max(p.cutoffs)>=1 | min(p.cutoffs)<=0 | is.unsorted(p.cutoffs))
	{
		return(cat("Error: The values of p.cutoffs must be greater than 0, less than 1, and in ascending order.","\n"))
	}	
	if(length(pal)==1 & is.character(pal))
	{
		if(!is.element(pal, rownames(brewer.pal.info)))
		{
			return(cat("Error: pal must be either a palette from RColorBrewer(to see a list: rownames(brewer.pal.info)) or a vector of colors.","\n"))
		}
	}else if (length(pal)!=length(p.cutoffs)+1)
	{
		return(cat("Error: The numbers of colors in pal must be one more than the number of values in p.cutoffs.","\n"))
	}
	if(test!="Stouffers" & test!="Fishers")
	{
		return(cat("Error: Value of test must be either \"Stouffers\" or \"Fishers\".","\n"))
	}
	#if(test=="Stouffers" & side==2)
	#{
    #		cat("Caution: Stouffer's Method is designed for 1-sided p-values.", "\n")
	#}
	if(test=="Fishers" & side==1)
	{
		cat("Caution: For Fisher's Method applied to one-tailed p-values, significance thresholds for small p-values (near 0) are more meaningful than for large p-values (near 1).", "\n\n")
	}
	if(branch.label!=TRUE & branch.label!=FALSE)
	{
		return(cat("Error: Value of branch.label must be either TRUE or FALSE.","\n"))
	}
	if(tip.color!=TRUE & tip.color!=FALSE)
	{
		return(cat("Error: Value of tip.color must be either TRUE or FALSE.","\n"))
	}
	if(edge.color!=TRUE & edge.color!=FALSE)
	{
		return(cat("Error: Value of edge.color must be either TRUE or FALSE.","\n"))
	}
	if(!is.numeric(tip.label.size))
	{
		return(cat("Error: Value of tip.label.size must be numeric"))
	}
	if(!is.numeric(branch.label.size))
	{
		return(cat("Error: Value of branch.label.size must be numeric"))
	}
	if(use.edge.length!=TRUE & use.edge.length!=FALSE)
	{
		return(cat("Error: Value of use.edge.length must be either TRUE or FALSE.","\n"))
	}
	if(branch != "node" & branch != "edge")
	{
		return(cat("Error: Value of branch must be either 'node' or 'edge'.","\n"))
	}
	if(root.edge!=TRUE & root.edge!=FALSE)
	{
		return(cat("Error: Value of root.edge must be either TRUE or FALSE.","\n"))
	}
	if(root.edge & branch=="node")
	{
	   cat("Warning: When root.edge=TRUE and branch='node', root edge color is not significance-based.","\n")
	}

	#Get tipcolor and edgecolor from tip.colors and edge.colors, respectively.  These
	#  are matrices of color values used in the plotphylo2 function
	tipcolor <- tip.colors(tree, unsorted.pvalues, p.cutoffs, pal, test, adjust, 
		side, method)
	edgecolor <- edge.colors(tree, unsorted.pvalues, p.cutoffs, pal, test, adjust, 
		side, method, branch)
	#Test whether the tip.color and edge.color parameters are true.  If tip.color (the
	#  parameter from our function) is true, then tip.color=tipcolor in plot (this is
	#  actually plotphylo2).  If false, then tip.color="black."  If edge.color (the
	#  parameter from our function) is true, then edge.color=edgecolor.  If false, then
	#  edge.color="black." We plot our tree.  We set show.node.label to false.  If the
	#  user wants node labels, we handle it by the nodelabels() function below.  type
	#  means the type of plot.  cex is the size of the tip labels and is set by the
	#  parameter tip.label.size.
	{if(tip.color & edge.color)	plotphylo2(tree, tip.color=tipcolor, show.node.label=FALSE, 
		edge.color=edgecolor, show.tip.label=TRUE, type=type, cex=tip.label.size,
		use.edge.length=use.edge.length, edge.width=edge.width, root.edge=root.edge)
	else if(tip.color & !edge.color)	plotphylo2(tree, tip.color=tipcolor, 
		show.node.label=FALSE, edge.color="black", show.tip.label=TRUE, type=type, 
		cex=tip.label.size, use.edge.length=use.edge.length, edge.width=edge.width, root.edge=root.edge)
	else if(!tip.color & edge.color)	plotphylo2(tree, tip.color="black", 
		show.node.label=FALSE, edge.color=edgecolor, show.tip.label=TRUE, 
		type=type, cex=tip.label.size, use.edge.length=use.edge.length, edge.width=edge.width, root.edge=root.edge)
	else	plotphylo2(tree, tip.color="black", show.node.label=FALSE, edge.color="black", 
		show.tip.label=TRUE, type=type, cex=tip.label.size,
		use.edge.length=use.edge.length, edge.width=edge.width, root.edge=root.edge)}
	#if the parameter branch.label is true, then we add labels to the branches.
	#  cex is the size of the branch labels and is set by the parameter branch.label.size.
	#  If we do not want to have frames around these labels, then set frame to "none." 
	#  Note that when branch='edge', the tip edges are not labeled (to avoid clutter).
	if(branch.label==TRUE){
	   br.frame <- match.arg(branch.label.frame, c("rect", "circle", "none"))
	   if(branch=='node'){ nodelabels(bg='white',cex=branch.label.size, frame=br.frame) }
	   if(branch=='edge'){ 	
			n.tips <- num.tips(tree)
			n.total.nodes <- num.total.nodes(tree)
			Branch <- as.character(c(tree$tip.label,(n.tips+1):n.total.nodes))
			Branch[n.tips+1] <- 'root'
			f1 <- data.frame(Node=Branch[-c(1:(n.tips+1))])
			f2 <- data.frame(Node=tree$edge[,2], Edge=c(1:num.edges(tree)))
			f <- merge(f1,f2, sort=F)
			edgelabels(edge=f$Edge,bg='white',cex=branch.label.size, frame=br.frame)
           }		
		}
}
