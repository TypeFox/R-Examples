
#Functions that accept the tree and return the number of edges=num.edges; number of
#  tips=num.tips
#number of internal nodes=num.internal.nodes; total number of nodes (internal + tip) =
#  num.nodes
#Internal functions
#Arguments: One or more of the functions plotSigTree, export.inherit, and
#	export.figtree call these functions (directly or indirectly).  Find
#	argument descriptions in plotSigTree, export.inherit, or export.figtree.
num.edges <- function(tree) length(tree$edge[,1])
num.tips <- function(tree) length(tree$tip.label)
num.internal.nodes <- function(tree) tree$Nnode
num.total.nodes <- function(tree) length(tree$tip.label)+tree$Nnode

srt.pvalues <- function(tree, unsorted.pvalues)
{
#The order of p value labels is likely different than the order of the tip labels in
#  the tree.
#This sorts the original p values according to the order of the tips in the tree using
#  the merge function.
#It returns the p values in sorted order (sorted in the same way as the order of the
#  tree tips).
#Internal function
#Arguments: One or more of the functions plotSigTree, export.inherit, and
#	export.figtree call this function (directly or indirectly).  Find
#	argument description in plotSigTree, export.inherit, or export.figtree.
	sorted.pvalues <- merge(tree$tip.label, unsorted.pvalues, by.x=1, by.y=1, all=TRUE, sort=FALSE) 
	return(sorted.pvalues)	
}

stouffers <- function(x)
{	
#Function that accepts a vector of values and calculates a family-wide p value based on
#  Stouffer's method
#Internal function
#Argument: x is a vector of p-values
	pnorm(sum(qnorm(x)) /sqrt(length(x)))
}

fishers <- function(y)
{
#Function that accepts a vector of values and calculates a family-wide p value based on
#  Fisher's method
#Internal function
#Argument: y is a vector of numbers	
	pchisq((-2)*(sum(log(y))),2*length(y),lower.tail=FALSE)
}

index.matrix <- function(tree)
{
#Create an index data frame called index.  We do this so we can know which descendants
#  belong to each node, or family (used later to calculate family wide pvalues later).
#  Each column corresponds to a node (family).  The rows correspond to the tips.  If
#  there is a 1 in a cell, then the row # (tip) belongs to the column # (node/family).
#  For example, if column 2, row 2 had value 1, then tip 2 (tree$tip.label[2] belongs
#  to node/family 2.  If it had value 0, it doesn't belong to node/family 2.
#This assumes that each internal nodes all have smaller numbers than each of their
#  descendants, not including tips.
#Internal function
#Argument: One or more of the functions plotSigTree, export.inherit, and
#	export.figtree call this function (directly or indirectly).  Find
#	argument description in plotSigTree, export.inherit, or export.figtree.

	#Call num.tips, num.internal.nodes, and num.total.nodes to get us the # of edges, #
	#  of internal nodes, and # of total nodes, respectively.
	n.tips <- num.tips(tree)
	n.internal.nodes <- num.internal.nodes(tree)
	n.total.nodes <- num.total.nodes(tree)
	n.edges <- num.edges(tree)
	#Create an identity matrix, tipmatrix, corresponding to the first (n.tips) columns.
	#  This is done because the first columns (columns 1 to n.tips) correspond only to
	#  tips and not internal nodes and will only include themselves in their family.
	#  For example, if there were 25 tips, then the index matrix would contain 25 rows
	#  and 25 columns.  Column 12 would only have a 1 in row 12, with the rest being
	#  0's.  After (n.tips), the nodes are internal nodes (not tips), and hence will
	#  have more than one tip in their family.  We create internalnodesmatrix, simply a
	#  matrix full of 0's with (n.tips) rows and (n.internal.nodes) columns.  Then
	#  append internalnodesmatrix to tipmatrix to get index, a (currently) blank matrix
	#  of size (n.tips) rows and (n.tips + n.internal.nodes) columns
	diag(n.tips) -> tipmatrix
	mat.or.vec(n.tips, n.internal.nodes) -> internalnodesmatrix
	index <- cbind(tipmatrix,internalnodesmatrix)
	#Sort the tree$edge matrix by the first column (in descending order) and call it
	#  tree$edgesort in order for the following for loop to capture all of what we need.
	#  *This assumes that our tree has its nodes (and edge matrix) assigned the way we
	#  expect.  This is that the root is the (n.tips+1) node and that an internal node
	#  always has a lower number than any of its descendants.  This is important so that
	#  our sort in the following line of code will work.  The edge matrix is a matrix
	#  that identifies how all of the nodes are connected to each other.  The row number
	#  corresponds to the edge number.  For examples, if in our edge matrix, row 2 has
	#  value 12 in column 1 and value 13 in column 2, then edge 2 connects node 12 to
	#  node 13.  For a given row, we assume that the edge matrix always has a higher
	#  number in column 2 if column 2 is a tip. If column 2 is a tip, we assume column
	#  2 will have a higher number.  We also assume that tips can only be in column 2.
	tree$edge[order(-tree$edge[,1]),] -> tree$edgesort
	#Fill in the index matrix with 1's if the tip corresponding to the row belongs to
	#  the family corresponding to the column.  The way tree$edgesort is sorted allows
	#  us to use the following for loop to do this. Our tree$edgesort matrix still holds
	#  the information about which nodes are connected together.  The edges are just
	#  sorted in a different way now. Because our index matrix has the identity matrix as
	#  the first (n.tips x n.tips) rows x columns, and because the identity matrix contains
	#  the information that the tips belong to their own family, and because the
	#  tree$edgesort matrix is sorted high to low, our for loop fills in the index matrix
	#  with 1's for all of the rows/tips that are descenants of a given column.  During the
	#  first iteration, the tree$edgesort[1,1]th (1) column of index will be added to the
	#  tree$edgesort[1,2]th (2) column of index and assigned to the tree$edgesort[1,1]th
	#  column of index.  Because (2) is a tip, it will have exactly one '1.'  Now, (1),
	#  which is an ancestor of (2), will have one '1.' During the second iteration, (1)
	#  will get another '1,' bringing its total to 2.  Thus (1) is an ancestor of both (2)
	#  and (3): the tree$edgesort[2,2]th (1) column of index.  This will continue to work
	#  because of the way it is sorted and the assumption that for non-tips, the ancestor
	#  has a smaller node number than the descendant.
#also assumes that each node only has 2 children
	for(i in 1:(n.edges)) #changed from n.total.nodes-1 to n.edges
	{			
		index[,tree$edgesort[i,1]] <- index[,tree$edgesort[i,1]] + index[,tree$edgesort[i,2]]
	}
	return(index)
}



p.p2.ADJ.p1 <- function(p, method)
{
#Function to convert a vector of one-tailed p-values to two-tailed, perform p-value
#  adjustment on the p-values (for multiple hypothesis testing), and then go back to
#  the one-tailed [adjusted] scale.
#This assumes null Mean2=Mean1 and alt: Mean2>Mean1
#P-value adjustment methods are found in ?p.adjust
#Argument: p is a vector of one-tailed p-values
#Arguments: One or more of the functions plotSigTree, export.inherit, and
#	export.figtree call this function (directly or indirectly).  Find
#	argument descriptions in plotSigTree, export.inherit, or export.figtree.
  	t <- p <= 0.5 # T/F of left tail
  	p1 <- p2 <- rep(NA,length(p))
  	p2[t] <- 2*p[t]
  	p2[!t] <- 2*(1-p[!t])
  	p.ADJ <- p.adjust(p2,method=method)
  	p1[t] <- p.ADJ[t]/2
  	p1[!t] <- 1-p.ADJ[!t]/2
  	return(p1)
}




result <- function(tree, unsorted.pvalues, test, adjust, side, method)
{
#Calculate pvalues based on Stouffer's or Fisher's tests.  Return results, a data frame
#  that contains this information.
#Internal function
#Arguments: One or more of the functions plotSigTree, export.inherit, and
#	export.figtree call this function (directly or indirectly).  Find
#	argument descriptions in plotSigTree, export.inherit, or export.figtree.

	#Call other functions to get n.total.nodes, the index matrix, sorted pvalues, and
	#  create an empty data frame called results
	n.total.nodes <- num.total.nodes(tree)
	index <- index.matrix(tree)
	results <- data.frame()	
	sorted.pvalues <- srt.pvalues(tree, unsorted.pvalues)
	#Fill in results data frame, which contains the p-values for all of the
	#  nodes/families for the selected test (Stouffer's or Fisher's).  temp contains the
	#  pvalues (obtained from the adjusted.pvalues data frame) for the jth column (node).
	#  temp contains only the p-values for the rows where index[,j]==1, ie. the tips
	#  belonging to the jth column.  We then call the Stouffer's and Fisher's functions
	#  on temp and assign the values we obtain to the first column of results.  Each row
	#  of results represents a different node.  
	if(test=="Stouffers")
	{
		for(j in 1:n.total.nodes)
		{
			temp <- sorted.pvalues[index[,j]==1,2]
			results[j,1] <- stouffers(temp)
			names(results)[1]<-"Stouffer's"
		}
	}else
	{
		for(j in 1:n.total.nodes)
		{
			temp <- sorted.pvalues[index[,j]==1,2]
			results[j,1] <- fishers(temp)
			names(results)[1]<-"Fisher's"
		}
	}
	#Now we need to apply the p-value adjustment if the adjust argument is TRUE.  If
	#	not, this step is skipped. Afterwards, results is returned.
	if(adjust==TRUE)
	{
		if(side==1)
		{
			results[,1] <- p.p2.ADJ.p1(results[,1], method)
		}else	
		{
			results[,1] <- p.adjust(results[,1], method=method)
		}
	}
	
	return(results)
}






tip.colors<-function(tree, unsorted.pvalues, p.cutoffs, pal, test, 
	adjust, side, method)
{
#Create tipcolor, which is a matrix of color values that determines how the tips/OTUs
#  are to be colored. It is of (n.tips) length and contains color values. Each row
#  represents a tip. 
#Internal function
#Arguments: One or more of the functions plotSigTree, export.inherit, and
#	export.figtree call this function (directly or indirectly).  Find
#	argument descriptions in plotSigTree, export.inherit, or export.figtree.

	#get results from the result function and create tipcolor, a matrix.
	results <- result(tree, unsorted.pvalues, test, adjust, side, method)
	matrix() -> tipcolor
	#create p.cutoffs.new, which is just p.cutoffs with 1 appended to the end.  n.cutoffs
	#  is the length of p.cutoffs.new.
	p.cutoffs.new<-c(p.cutoffs, 1)
	n.cutoffs <- length(p.cutoffs.new)
	#Test to see if pal is of length 1.  If it is, it will be a valid RColorBrewer
	#  palette due to the error checking in the function that called it (assuming that
	#  either plotSigTree or export.figtree called it).  If it is a valid palette, we
	#  create cols and assign it a vector of hexadecimal colors from RColorBrewer based
	#  on n.cutoffs and pal.  If it is not of length 1, then it is a vector of (we
	#  assume) valid colors.  Then we simply assign pal to cols.
	if(length(pal)==1)
	{	
			cols <- brewer.pal(n.cutoffs, pal)
			if(side==1){ cols[ceiling(n.cutoffs/2)] <- brewer.pal(7, "Greys")[2] }
	}else 
	{	
		cols <- pal
	}
	#n.tips is the number of tips in tree.
	n.tips <- num.tips(tree)
	#for loop that works its way backwards from n.cutoffs to 1. for each k, results[,1]
	#  <= p.cutoffs.new[k] is a vector of TRUE/FALSE values. For all of the TRUE values,
	#  cols[k] is assigned to tipcolor.  cols is a vector of colors.  For example, the
	#  first iteration is when k=n.cutoffs.  This is the length of p.cutoffs.new.  If our
	#  cutoffs were (.01, .05, .10, .90, .95, .99), then p.cutoffs.new would be (.01, .05,
	#  .10, .90, .95, .99, 1) and n.cutoffs would be 7.  Then k=7.  Cols would have 7
	#  colors in it.  cols[7], the 7th color, would be assigned to all of the rows of
	#  tipcolor where the corresponding p-value was less than or equal to 1.  As a matter
	#  of fact, all rows would be assigned this color because all of the p-values are by
	#  definition less than or equal to 1.  The second iteration would assign cols[6] to
	#  all of the rows where the p-values were <=.99 and so on.  The last iteration would
	#  only assign cols[1] to values that were <=.01.  Note that with the way we did our
	#  logic in this for loop, we needed to append 1 to p.cutoffs (and create
	#  p.cutoffs.new).  If we did not, we would have not assigned values to our top
	#  interval, for example (.99,1].
	for(k in n.cutoffs:1)
	{
		tipcolor[results[,1] <= p.cutoffs.new[k]]<-cols[k]
	}
	#Take only the first (n.tips) values of tipcolor and assign to tipcolor.  We do this
	#  because tipcolor is originally too long.  It should only be (n.tips) long because
	#  there are only (n.tips) tips.
	tipcolor <- tipcolor[1:n.tips]
	return(tipcolor)
}




edge.colors <- function(tree, unsorted.pvalues, p.cutoffs, pal, test, adjust, side, method, branch)
{
# This is an updated version of the edge.color function.
# It allows for root edge color calculation (at end of returned vector) -- for use when branch='edge'
#Create edgecolor, which is a matrix that determines how the edges/families are to be
#  colored.  It is of length (n.edges+1)-i think and contains color values.  Each row
#  represents an edge.
#Internal function
#Arguments: One or more of the functions plotSigTree, export.inherit, and
#	export.figtree call this function (directly or indirectly).  Find
#	argument descriptions in plotSigTree, export.inherit, or export.figtree.

	#get results from result function; create edgecolor, a matrix
	results <- result(tree, unsorted.pvalues, test, adjust, side, method)
	matrix() -> edgecolor
	#create p.cutoffs.new, which is just p.cutoffs with 1 appended to the end.  n.cutoffs
	#  is the length of p.cutoffs.new. cols is a vector of hexadecimal colors from
	#  RColorBrewer based on n.cutoffs and pal
	p.cutoffs.new <- c(p.cutoffs, 1)
	n.cutoffs <- length(p.cutoffs.new)
	
		#Test to see if pal is of length 1.  If it is, it will be a valid RColorBrewer
		#  palette due to the error checking in the function that called it (assuming that
		#  either plotSigTree or export.figtree called it).  If it is a valid palette, we
		#  create cols and assign it a vector of hexadecimal colors from RColorBrewer
		#  based on n.cutoffs and pal.  If it is not of length 1, then it is a vector of
		#  (we assume) valid colors.  Then we simply assign pal to cols.
	if(length(pal)==1)
	{	
			cols <- brewer.pal(n.cutoffs, pal)
			if(side==1){ cols[ceiling(n.cutoffs/2)] <- brewer.pal(7, "Greys")[2] }
	}else
	{
		cols <- pal
	}
		#For loop that works its way backwards from n.cutoffs to 1.
		#  (1) results[tree$edge[,1],1]<=p.cutoffs.new[k] is a vector of TRUE/FALSE
		#  values.  We use results[tree$edge[,1],1] in our comparison because we want to
		#  find the p-value corresponding to the left node (if looking at the tree where
		#  the root is on the left and the tips on the right; otherwise the more interior
		#  node; the ancestor) of the edge to determine the coloring for that edge.  This
		#  is why we choose tree$edge[,1] instead of tree$edge[,2].  (1) will be in the
		#  order of the edge matrix.  Thus edgecolor will be in the correct order.  
		#For example, on the first iteration, k=n.cutoffs.  This is the length of
		#  p.cutoffs.new.  If our cutoffs were (.01, .05, .10, .90, .95, .99), then
		#  p.cutoffs.new would be (.01, .05, .10, .90, .95, .99, 1) and n.cutoffs would be
		#  7.  Then k=7.  Cols would have 7 colors in it.  cols[7], the 7th color, would be
		#  assigned to all of the rows (edges) of edgecolor where the corresponding p-value
		#  (from results[tree$edge[,1],1]) was less than or equal to 1.  As a matter of
		#  fact, all rows would be assigned this color because all of the p-values are by
		#  definition less than or equal to 1.  The second iteration would assign cols[6]
		#  to all of the rows where the corresponding p-value were <=.99.  The last
		#  iteration would only assign cols[1] to values that were <=.01.  Note that with
		#  the way we did our logic in this for loop, we needed to append 1 to p.cutoffs
		#  (and create p.cutoffs.new).  If we did not, we would have not assigned values to
		#  our top interval, for example (.99,1].
	
    # Depending on definition of branch (node or edge), get edge colors
	# (Prior to package version 1.2, only "node" option was done.)
	if(branch=="node")
	  {
    	for(k in n.cutoffs:1)
	     {
		     edgecolor[results[tree$edge[,1],1]<=p.cutoffs.new[k]] <- cols[k]
	     }
	  }
	if(branch=="edge")
	  {
    	for(k in n.cutoffs:1)
	     {
		     edgecolor[results[tree$edge[,2],1]<=p.cutoffs.new[k]] <- cols[k]
	     }
		# Add final color for root edge (not included in tree$edge)
		ntips <- num.tips(tree)
		root.pval <- results[(ntips+1),1]
		t <- root.pval <= p.cutoffs.new
		edgecolor[nrow(results)] <- cols[t][1]
      }
	
	return(edgecolor)
}




# Function based on ape's plot.phylo function
# - all arguments are the same, except for root.edge.col
# - This function is for use in the SigTree package,
#   called by the plotSigTree function, allowing full branch coloring
#   (including the 'perpendicular-to-the-root' half-edges that were impossible to color
#   differently in plot.phylo)
# - the argument root.edge.col is a color to use for the root edge.
#   If root.edge.col is NULL and root.edge is TRUE, then the root edge
#   color will be either black (default) or else (if the length of
#   edge.color is one more than the number of [non-root] edges) it will be taken
#   as the last color in edge.color.
# - This function also allows for root.edge=TRUE for both type='phylogram'
#   and type='fan' (plot.phylo only supported type='phylogram')
# Beginning in SigTree version 1.5, the sections here relating to .C() calls
#   are updated to coincide with the newer version of SigTree.c (which itself
#   is a copy of the updated plot_phylo.c file from the ape package)
plotphylo2 <-
function (x, type = "phylogram", use.edge.length = TRUE, node.pos = NULL, 
    show.tip.label = TRUE, show.node.label = FALSE, edge.color = "black", 
    edge.width = 1, edge.lty = 1, font = 3, cex = par("cex"), 
    adj = NULL, srt = 0, no.margin = FALSE, root.edge = FALSE, 
    label.offset = 0, underscore = FALSE, x.lim = NULL, y.lim = NULL, 
    direction = "rightwards", lab4ut = "horizontal", tip.color = "black", 
    plot = TRUE, rotate.tree = 0, open.angle = 0, root.edge.col=NULL, ...) 
{
	if(root.edge && is.null(x$root.edge)){x <- compute.brlen(x); x$root.edge <- median(x$edge.length)}
    Ntip <- length(x$tip.label)
    if (Ntip < 2) {
        warning("found less than 2 tips in the tree")
        return(NULL)
    }
    if (any(tabulate(x$edge[, 1]) == 1)) 
        stop("there are single (non-splitting) nodes in your tree; you may need to use collapse.singles()")
    .nodeHeight <- function(Ntip, Nnode, edge, Nedge, yy)
        .C(node_height, as.integer(Ntip), as.integer(Nnode),
           as.integer(edge[, 1]), as.integer(edge[, 2]),
           as.integer(Nedge), as.double(yy), PACKAGE = "SigTree")[[6]]

    .nodeDepth <- function(Ntip, Nnode, edge, Nedge, node.depth)
        .C(node_depth, as.integer(Ntip), as.integer(Nnode),
           as.integer(edge[, 1]), as.integer(edge[, 2]),
           as.integer(Nedge), double(Ntip + Nnode), as.integer(node.depth), PACKAGE = "SigTree")[[6]]

    .nodeDepthEdgelength <- function(Ntip, Nnode, edge, Nedge, edge.length)
        .C(node_depth_edgelength, as.integer(Ntip),
           as.integer(Nnode), as.integer(edge[, 1]),
           as.integer(edge[, 2]), as.integer(Nedge),
           as.double(edge.length), double(Ntip + Nnode), PACKAGE = "SigTree")[[7]]
    Nedge <- dim(x$edge)[1]
    Nnode <- x$Nnode
    ROOT <- Ntip + 1
    type <- match.arg(type, c("phylogram", "cladogram", "fan", 
        "unrooted", "radial"))
    direction <- match.arg(direction, c("rightwards", "leftwards", 
        "upwards", "downwards"))
    if (is.null(x$edge.length)) 
        use.edge.length <- FALSE
    if (type %in% c("unrooted", "radial") || !use.edge.length || 
        is.null(x$root.edge) || !x$root.edge) 
        root.edge <- FALSE
    if (type == "fan" && root.edge) { use.edge.length <- FALSE }
    # Allow for custom coloring of root edge, using 'extra' color returned by edge.colors function when root.edge=TRUE
	if(root.edge & is.null(root.edge.col))
	  {
	    if(length(edge.color)==(Nedge+1)) {
		    root.edge.col <- rev(edge.color)[1] 
			 } else {root.edge.col <- 'black'}
	   }
	if(root.edge & is.null(root.edge.col) & length(edge.color)==(Nedge+1)){root.edge.col <- rev(edge.color)[1]}
    phyloORclado <- type %in% c("phylogram", "cladogram")
    horizontal <- direction %in% c("rightwards", "leftwards")
    xe <- x$edge
    if (phyloORclado) {
        phyOrder <- attr(x, "order")
        if (is.null(phyOrder) || phyOrder != "cladewise") {
            x <- reorder(x)
            if (!identical(x$edge, xe)) {
                ereorder <- match(x$edge[, 2], xe[, 2])
                if (length(edge.color) > 1) {
                  edge.color <- rep(edge.color, length.out = Nedge)
                  edge.color <- edge.color[ereorder]
                }
                if (length(edge.width) > 1) {
                  edge.width <- rep(edge.width, length.out = Nedge)
                  edge.width <- edge.width[ereorder]
                }
                if (length(edge.lty) > 1) {
                  edge.lty <- rep(edge.lty, length.out = Nedge)
                  edge.lty <- edge.lty[ereorder]
                }
            }
        }
        yy <- numeric(Ntip + Nnode)
        TIPS <- x$edge[x$edge[, 2] <= Ntip, 2]
        yy[TIPS] <- 1:Ntip
    }
    z <- reorder(x, order = "pruningwise")
    if (phyloORclado) {
        if (is.null(node.pos)) {
            node.pos <- 1
            if (type == "cladogram" && !use.edge.length) 
                node.pos <- 2
        }
        if (node.pos == 1) 
           { yy <- .nodeHeight(Ntip, Nnode, z$edge, Nedge, yy) }  else {
            ans <- .C(node_height_clado, as.integer(Ntip),
                    as.integer(Nnode), as.integer(z$edge[, 1]),
                    as.integer(z$edge[, 2]), as.integer(Nedge),
                    double(Ntip + Nnode), as.double(yy), PACKAGE = "SigTree")
            xx <- ans[[6]] - 1
            yy <- ans[[7]]
        }
        if (!use.edge.length) {
            if (node.pos != 2) 
                xx <- .nodeDepth(Ntip, Nnode, z$edge, Nedge) - 
                  1
            xx <- max(xx) - xx
        }  else {
            xx <- .nodeDepthEdgelength(Ntip, Nnode, z$edge, Nedge, 
                z$edge.length)
        }
    }  else {
        twopi <- 2 * pi
        rotate.tree <- twopi * rotate.tree/360
        switch(type, fan = {
            TIPS <- x$edge[which(x$edge[, 2] <= Ntip), 2]
            xx <- seq(0, twopi * (1 - 1/Ntip) - twopi * open.angle/360, 
                length.out = Ntip)
            theta <- double(Ntip)
            theta[TIPS] <- xx
            theta <- c(theta, numeric(Nnode))
            theta <- .nodeHeight(Ntip, Nnode, z$edge, Nedge, 
                theta)
            if (use.edge.length) {
                r <- .nodeDepthEdgelength(Ntip, Nnode, z$edge, 
                  Nedge, z$edge.length)
            } else {
                r <- .nodeDepth(Ntip, Nnode, z$edge, Nedge)
                r <- 1/r
            }
            theta <- theta + rotate.tree
            xx <- r * cos(theta)
            yy <- r * sin(theta)
        }, unrooted = {
            nb.sp <- .nodeDepth(Ntip, Nnode, z$edge, Nedge)
            XY <- if (use.edge.length) unrooted.xy(Ntip, Nnode, 
                z$edge, z$edge.length, nb.sp, rotate.tree) else unrooted.xy(Ntip, 
                Nnode, z$edge, rep(1, Nedge), nb.sp, rotate.tree)
            xx <- XY$M[, 1] - min(XY$M[, 1])
            yy <- XY$M[, 2] - min(XY$M[, 2])
        }, radial = {
            X <- .nodeDepth(Ntip, Nnode, z$edge, Nedge)
            X[X == 1] <- 0
            X <- 1 - X/Ntip
            yy <- c((1:Ntip) * twopi/Ntip, rep(0, Nnode))
            Y <- .nodeHeight(Ntip, Nnode, z$edge, Nedge, yy)
            xx <- X * cos(Y + rotate.tree)
            yy <- X * sin(Y + rotate.tree)
        })
    }
    if (phyloORclado) {
        if (!horizontal) {
            tmp <- yy
            yy <- xx
            xx <- tmp - min(tmp) + 1
        }
        if (root.edge) {
            if (direction == "rightwards") 
                xx <- xx + x$root.edge
            if (direction == "upwards") 
                yy <- yy + x$root.edge
        }
    }
    if (no.margin) 
        par(mai = rep(0, 4))
    if (is.null(x.lim)) {
        if (phyloORclado) {
            if (horizontal) {
                x.lim <- c(0, NA)
                pin1 <- par("pin")[1]
                strWi <- strwidth(x$tip.label, "inches")
                xx.tips <- xx[1:Ntip] * 1.04
                alp <- try(uniroot(function(a) max(a * xx.tips + 
                  strWi) - pin1, c(0, 1e+06))$root, silent = TRUE)
                if (is.character(alp)) 
                  { tmp <- max(xx.tips) * 1.5 } else {
                  tmp <- if (show.tip.label) 
                   { max(xx.tips + strWi/alp) }  else max(xx.tips)
                }
                x.lim[2] <- tmp
            } else x.lim <- c(1, Ntip)
        } else switch(type, fan = {
            if (show.tip.label) {
                offset <- max(nchar(x$tip.label) * 0.018 * max(yy) * 
                  cex)
                x.lim <- c(min(xx) - offset, max(xx) + offset)
            } else x.lim <- c(min(xx), max(xx))
        }, unrooted = {
            if (show.tip.label) {
                offset <- max(nchar(x$tip.label) * 0.018 * max(yy) * 
                  cex)
                x.lim <- c(0 - offset, max(xx) + offset)
            } else x.lim <- c(0, max(xx))
        }, radial = {
            if (show.tip.label) {
                offset <- max(nchar(x$tip.label) * 0.03 * cex)
                x.lim <- c(-1 - offset, 1 + offset)
            } else x.lim <- c(-1, 1)
        })
    } else if (length(x.lim) == 1) {
        x.lim <- c(0, x.lim)
        if (phyloORclado && !horizontal) 
            x.lim[1] <- 1
        if (type %in% c("fan", "unrooted") && show.tip.label) 
            x.lim[1] <- -max(nchar(x$tip.label) * 0.018 * max(yy) * 
                cex)
        if (type == "radial") 
            x.lim[1] <- if (show.tip.label) 
                { -1 - max(nchar(x$tip.label) * 0.03 * cex) } else { -1 }
    }
    if (phyloORclado && direction == "leftwards") 
        xx <- x.lim[2] - xx
    if (is.null(y.lim)) {
        if (phyloORclado) {
            if (horizontal) 
               { y.lim <- c(1, Ntip) }  else {
                y.lim <- c(0, NA)
                pin2 <- par("pin")[2]
                strWi <- strwidth(x$tip.label, "inches")
                yy.tips <- yy[1:Ntip] * 1.04
                alp <- try(uniroot(function(a) max(a * yy.tips + 
                  strWi) - pin2, c(0, 1e+06))$root, silent = TRUE)
                if (is.character(alp)) 
                 { tmp <- max(yy.tips) * 1.5 }  else {
                  tmp <- if (show.tip.label) 
                    { max(yy.tips + strWi/alp)} else max(yy.tips)
                }
                y.lim[2] <- tmp
            }
        } else switch(type, fan = {
            if (show.tip.label) {
                offset <- max(nchar(x$tip.label) * 0.018 * max(yy) * 
                  cex)
                y.lim <- c(min(yy) - offset, max(yy) + offset)
            } else y.lim <- c(min(yy), max(yy))
        }, unrooted = {
            if (show.tip.label) {
                offset <- max(nchar(x$tip.label) * 0.018 * max(yy) * 
                  cex)
                y.lim <- c(0 - offset, max(yy) + offset)
            } else y.lim <- c(0, max(yy))
        }, radial = {
            if (show.tip.label) {
                offset <- max(nchar(x$tip.label) * 0.03 * cex)
                y.lim <- c(-1 - offset, 1 + offset)
            } else y.lim <- c(-1, 1)
        })
    } else if (length(y.lim) == 1) {
        y.lim <- c(0, y.lim)
        if (phyloORclado && horizontal) 
            y.lim[1] <- 1
        if (type %in% c("fan", "unrooted") && show.tip.label) 
            y.lim[1] <- -max(nchar(x$tip.label) * 0.018 * max(yy) * 
                cex)
        if (type == "radial") 
            y.lim[1] <- if (show.tip.label) 
               { -1 - max(nchar(x$tip.label) * 0.018 * max(yy) * 
                  cex) } else  -1
    }
    if (phyloORclado && direction == "downwards") 
        yy <- max(yy) - yy
    if (phyloORclado && root.edge) {
        if (direction == "leftwards") 
            x.lim[2] <- x.lim[2] + x$root.edge
        if (direction == "downwards") 
            y.lim[2] <- y.lim[2] + x$root.edge
    }
    asp <- if (type %in% c("fan", "radial", "unrooted")) 
        { 1 } else NA
    plot(0, type = "n", xlim = x.lim, ylim = y.lim, ann = FALSE, 
        axes = FALSE, asp = asp, ...)
	if (plot) {
      if (is.null(adj)) 
            adj <- if (phyloORclado && direction == "leftwards") 
               { 1 } else 0
        if (phyloORclado && show.tip.label) {
            MAXSTRING <- max(strwidth(x$tip.label, cex = cex))
            loy <- 0
            if (direction == "rightwards") {
                lox <- label.offset + MAXSTRING * 1.05 * adj
            }
            if (direction == "leftwards") {
                lox <- -label.offset - MAXSTRING * 1.05 * (1 - 
                  adj)
            }
            if (!horizontal) {
                psr <- par("usr")
                MAXSTRING <- MAXSTRING * 1.09 * (psr[4] - psr[3])/(psr[2] - 
                  psr[1])
                loy <- label.offset + MAXSTRING * 1.05 * adj
                lox <- 0
                srt <- 90 + srt
                if (direction == "downwards") {
                  loy <- -loy
                  srt <- 180 + srt
                }
            }
        }
		
        if (type == "phylogram") {
            phylogram.plot(x$edge, Ntip, Nnode, xx, yy, horizontal, 
                edge.color, edge.width, edge.lty)
        } else {
            if (type == "fan") {

                ereorder <- match(z$edge[, 2], x$edge[, 2])
                if (length(edge.color) > 1) {
                  edge.color <- rep(edge.color, length.out = Nedge)
                  edge.color <- edge.color[ereorder]
                }
                if (length(edge.width) > 1) {
                  edge.width <- rep(edge.width, length.out = Nedge)
                  edge.width <- edge.width[ereorder]
                }
                if (length(edge.lty) > 1) {
                  edge.lty <- rep(edge.lty, length.out = Nedge)
                  edge.lty <- edge.lty[ereorder]
                }
                circularplot2(z$edge, Ntip, Nnode, xx, yy, theta, 
                  r, edge.color, edge.width, edge.lty, root.edge, root.edge.col)
				root.edge <- FALSE
            } else cladogram.plot(x$edge, xx, yy, edge.color, edge.width, 
                edge.lty)
        }
        if (root.edge){ 
            switch(direction, 
			       rightwards = segments(0, yy[ROOT], x$root.edge, yy[ROOT], col=root.edge.col, lwd=edge.width), 
				   leftwards = segments(xx[ROOT], yy[ROOT], xx[ROOT] + x$root.edge, yy[ROOT], col=root.edge.col, lwd=edge.width), 
                   upwards = segments(xx[ROOT], 0, xx[ROOT], x$root.edge, col=root.edge.col, lwd=edge.width), 
                   downwards = segments(xx[ROOT], yy[ROOT], xx[ROOT], yy[ROOT] + x$root.edge, col=root.edge.col, lwd=edge.width)
				   ) }
        if (show.tip.label) {
            if (is.expression(x$tip.label)) 
                underscore <- TRUE
            if (!underscore) 
                x$tip.label <- gsub("_", " ", x$tip.label)
            if (phyloORclado) 
                text(xx[1:Ntip] + lox, yy[1:Ntip] + loy, x$tip.label, 
                  adj = adj, font = font, srt = srt, cex = cex, 
                  col = tip.color)
            if (type == "unrooted") {
                if (lab4ut == "horizontal") {
                  y.adj <- x.adj <- numeric(Ntip)
                  sel <- abs(XY$axe) > 0.75 * pi
                  x.adj[sel] <- -strwidth(x$tip.label)[sel] * 
                    1.05
                  sel <- abs(XY$axe) > pi/4 & abs(XY$axe) < 0.75 * 
                    pi
                  x.adj[sel] <- -strwidth(x$tip.label)[sel] * 
                    (2 * abs(XY$axe)[sel]/pi - 0.5)
                  sel <- XY$axe > pi/4 & XY$axe < 0.75 * pi
                  y.adj[sel] <- strheight(x$tip.label)[sel]/2
                  sel <- XY$axe < -pi/4 & XY$axe > -0.75 * pi
                  y.adj[sel] <- -strheight(x$tip.label)[sel] * 
                    0.75
                  text(xx[1:Ntip] + x.adj * cex, yy[1:Ntip] + 
                    y.adj * cex, x$tip.label, adj = c(adj, 0), 
                    font = font, srt = srt, cex = cex, col = tip.color)
                } else {
                  adj <- abs(XY$axe) > pi/2
                  srt <- 180 * XY$axe/pi
                  srt[adj] <- srt[adj] - 180
                  adj <- as.numeric(adj)
                  xx.tips <- xx[1:Ntip]
                  yy.tips <- yy[1:Ntip]
                  if (label.offset) {
                    xx.tips <- xx.tips + label.offset * cos(XY$axe)
                    yy.tips <- yy.tips + label.offset * sin(XY$axe)
                  }
                  font <- rep(font, length.out = Ntip)
                  tip.color <- rep(tip.color, length.out = Ntip)
                  cex <- rep(cex, length.out = Ntip)
                  for (i in 1:Ntip) text(xx.tips[i], yy.tips[i], 
                    cex = cex[i], x$tip.label[i], adj = adj[i], 
                    font = font[i], srt = srt[i], col = tip.color[i])
                }
            }
            if (type %in% c("fan", "radial")) {
                xx.tips <- xx[1:Ntip]
                yy.tips <- yy[1:Ntip]
                angle <- atan2(yy.tips, xx.tips)
                if (label.offset) {
                  xx.tips <- xx.tips + label.offset * cos(angle)
                  yy.tips <- yy.tips + label.offset * sin(angle)
                }
                s <- xx.tips < 0
                angle <- angle * 180/pi
                angle[s] <- angle[s] + 180
                adj <- as.numeric(s)
                font <- rep(font, length.out = Ntip)
                tip.color <- rep(tip.color, length.out = Ntip)
                cex <- rep(cex, length.out = Ntip)
                for (i in 1:Ntip) text(xx.tips[i], yy.tips[i], 
                  x$tip.label[i], font = font[i], cex = cex[i], 
                  srt = angle[i], adj = adj[i], col = tip.color[i])
            }
        }
        if (show.node.label){ 
            text(xx[ROOT:length(xx)] + label.offset, yy[ROOT:length(yy)], 
                x$node.label, adj = adj, font = font, srt = srt, 
                cex = cex)
			}
    }
    L <- list(type = type, use.edge.length = use.edge.length, 
        node.pos = node.pos, show.tip.label = show.tip.label, 
        show.node.label = show.node.label, font = font, cex = cex, 
        adj = adj, srt = srt, no.margin = no.margin, label.offset = label.offset, 
        x.lim = x.lim, y.lim = y.lim, direction = direction, 
        tip.color = tip.color, Ntip = Ntip, Nnode = Nnode)
    assign("last_plot.phylo", c(L, list(edge = xe, xx = xx, yy = yy)), 
        envir = .PlotPhyloEnv)
    invisible(L)
}


## This modifies ape's circular.plot function so 'perpendicular-to-the-root' edges
## aren't all black when edge.color is a vector in type='fan' plot.
## This is called by SigTree's plotphylo2 when type='fan'
circularplot2 <- 
function (edge, Ntip, Nnode, xx, yy, theta, r, edge.color, edge.width, 
    edge.lty, root.edge, root.edge.col) 
{
    r0 <- r[edge[, 1]]
    r1 <- r[edge[, 2]]
    theta0 <- theta[edge[, 2]]
    costheta0 <- cos(theta0)
    sintheta0 <- sin(theta0)
    x0 <- r0 * costheta0
    y0 <- r0 * sintheta0
    x1 <- r1 * costheta0
    y1 <- r1 * sintheta0

	segments(x0, y0, x1, y1, col = edge.color, lwd = edge.width, 
        lty = edge.lty)
		
	tmp <- which(diff(edge[, 1]) != 0)
    start <- c(1, tmp + 1)
    Nedge <- dim(edge)[1]
    end <- c(tmp, Nedge)
	
 	   ## The 'co' object returned by foo function in circular.plot function
	   ## is a vector of colors whose length is 'Nnode' 
	   ## (the number of interior nodes, or the number of 'perpendicular-to-the-root' edges).
	   ## I need it to instead be of length 2*Nnode, with colors
       ## in appropriate order to color each half of the 'perpendicular-to-the-root' edges.
       ## So I'll need a 'foo2' function and also I'll need 
       ## the 'for' loop below to do half-edges.

    foo2 <- function(edge.feat, default){
        if (length(edge.feat) == 1){ 
            return(rep(edge.feat, 2*Nnode)) }  else {
            edge.feat <- rep(edge.feat, length.out = Nedge)
             # edge.feat should now be a vector of colors (or features)
             # in the same order as the edges of the tree;
             # edge.feat[k] is the color (feature) of edge k, and
             # so should be the color (feature) of the half of the
             # 'perpendicular-to-the-root' edge touching that edge (and 
             # heading towards the root) -- so make this
			 # be the returned feat.arc object
			feat.arc <- edge.feat
		   }
        feat.arc
    }
	co2 <- foo2(edge.color, "red")
	lw2 <- foo2(edge.width,1)
	ly2 <- foo2(edge.lty,1)

    for (k in 1:(Nnode)) {
        i <- start[k]
        j <- end[k]
		# circular.plot function defined:
		#    node k lies on 'perpendicular-to-the-root' edge joining edges i and j (i<j)
		# Here, we want to do this by halves of the 'perpendicular-to-the-root' edge.
		ngrad <- 100
        X <- rep(r[edge[i, 1]], ngrad)
        Y <- seq(theta[edge[i, 2]], theta[edge[j, 2]], length.out = ngrad)
		X1 <- X[1:(ngrad/2)]
		Y1 <- Y[1:(ngrad/2)]
		X2 <- X[(ngrad/2+1):ngrad]
		Y2 <- Y[(ngrad/2+1):ngrad]
        lines(X1 * cos(Y1), X1 * sin(Y1), col = co2[2*k-1], lwd = lw2[2*k-1], lty = ly2[2*k-1])
        lines(X2 * cos(Y2), X2 * sin(Y2), col = co2[2*k], lwd = lw2[2*k], lty = ly2[2*k])

		}

	# Now - artificially add root edge by connecting midpoint of that last 'perpendicular-to-the-root' edge
	# to the origin
	if(root.edge==TRUE){ 
	    lines(c(X2[1]*cos(Y2[1]), 0), c(X2[1]*sin(Y2[1]), 0), col=root.edge.col, lwd=lw2[Nnode], lty=ly2[Nnode])
	 }
		
}


