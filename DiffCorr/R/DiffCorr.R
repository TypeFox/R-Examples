##' Analyze and visualize differential correlations in biological networks
##'
##' \tabular{ll}{
##' Package: \tab DiffCorr \cr
##' Type: \tab Package \cr
##' Version: \tab 0.4.1 \cr
##' Date: \tab 2015-03-31 \cr
##' Depends: \tab igraph, pcaMethods, fdrtool \cr
##' License: \tab GPL (>=3) \cr
##' LazyLoad: \tab yes \cr
##' }
##'
##''
##' @name DiffCorr
##' @aliases DiffCorr
##' @docType package
##' @import fdrtool
##' @import igraph
##' @import multtest
##' @import pcaMethods
##' @title Differential correlations in omics datasets
##' @author Atsushi Fukushima, Kozo Nishida
NULL



##' The pre-treatment methods
##' 
##' @title scalingMethods
##' @param data a data matrix ([data.frame object] row: molecules,
##' col: samples or replicates)
##' @param methods the chosen methods.
##' @return the resulting data frame (or scaled data matrix)
##' @examples
##' scalingMethods(iris[,1:4], "level")
##' @author Atsushi Fukushima
##' @export
scalingMethods <- function(data, 
            methods = c("auto", "range", "pareto", "vast", "level", "power")) {
  methods <- match.arg(methods)
  if(ncol(data) > 1) {
    switch(methods,
      auto = {
        res <- apply(data, 1,
                         function(x) (x-mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE))
        return(data.frame(t(res), check.names=FALSE))
      },
      range = {
        res <- apply(data, 1,
                         function(x)
                         (x - mean(x, na.rm=TRUE))/(range(x)[2]-range(x)[1]))
        return(data.frame(t(res), check.names=FALSE))
      },
      pareto = {
        res <- apply(data, 1,
                         function(x)
                         (x - mean(x, na.rm=TRUE))/sqrt(sd(x, na.rm=TRUE)))
        return(data.frame(t(res), check.names=FALSE))
      },
      vast = {
        res <- apply(data, 1,
                         function(x) mean(x, na.rm=TRUE) *
                         (x - mean(x, na.rm=TRUE)) / (sd(x, na.rm=TRUE)**2))
        return(data.frame(t(res), check.names=FALSE))
      },
      level = {
        res <- apply(data, 1,
                         function(x) (x - mean(x, na.rm=TRUE)) / mean(x, na.rm=TRUE))
        return(data.frame(t(res), check.names=FALSE))
      },
      power = {
        res <- apply(data, 1,
                         function(x) sqrt(x) - mean(sqrt(x)))
        return(data.frame(t(res), check.names=FALSE))
      })
  }
}

##' Correlation Test
##'
##' @title Correlation Test
##' @param n the number of samples
##' @param r the correlation coefficient
##' @param method "pearson" and "spearman" can be used.
##' @return p-value
##' @examples
##' cor2.test(30, 0.6)
##' @references http://aoki2.si.gunma-u.ac.jp/R/cor2.html
##' @author Atsushi Fukushima
##' @export
cor2.test <- function(n, r, method = c("pearson", "kendall", "spearman")) {
  method <- match.arg(method)
  if (method == "pearson" || method == "spearman") {
    t <- abs(r) * sqrt( (n - 2) / (1 - r**2) )
    df <- n - 2
    p <- pt(t, df, lower.tail = FALSE) * 2
    c(p)
  }
  else {
    z <- abs(r) / sqrt( (4 * n + 10) / (9 * n * (n-1) ) )
    p <- pnorm(z, lower.tail = FALSE) * 2
    c(p)
  }
}



##' Compare two correlation coefficients using Fisher's Z-transformation
##' 
##' @title Compare two correlation coefficients
##' @param n1 sample size under condition 1
##' @param r1 correlation coefficient under condition 1
##' @param n2 sample size under condition 2
##' @param r2 correlation coefficient under condition 1
##' @return list of result (diff and p-value)
##' @examples
##' compcorr(10, 0.1, 10, 0.9)
##' @references
##' http://www.fon.hum.uva.nl/Service/Statistics/Two_Correlations.html
##' http://support.sas.com/ctx/samples/index.jsp?sid=494
##' http://support.sas.com/ctx/samples/index.jsp?sid=494
##' @author Atsushi Fukushima
##' @export
compcorr <- function(n1, r1, n2, r2){
  # Fisher's Z-transformation
  # ad hoc process
  num1a <- which(r1 >= 0.99)
  num2a <- which(r2 >= 0.99)
  r1[num1a] <- 0.99
  r2[num2a] <- 0.99
  num1b <- which(r1 <= -0.99)
  num2b <- which(r2 <= -0.99)
  r1[num1b] <- -0.99
  r2[num2b] <- -0.99
  # z
  z1 <- atanh(r1)
  z2 <- atanh(r2)
  # difference
  dz <- (z1 - z2)/sqrt(1/(n1 - 3) + (1/(n2 - 3)))
  # p-value
  pv <- 2*( 1 - pnorm(abs(dz)) )
  return(list(diff=dz, pval=pv))
}

##' Getting local false discovery rate (lfdr) using 'fdrtool' package
##' 
##' @title Getting local false discovery rate (lfdr)
##' @param r a vector of correlation coefficient under condition
##' @return list of lfdr
##' @examples
##' library("fdrtool")
##' data(pvalues)
##' get.lfdr(pvalues)
##' @references
##' Strimmer, K. Bioinformatics (2008) 24, 1461-1462
##' @author Atsushi Fukushima
##' @export
get.lfdr <- function(r) {
  fdr.out <- fdrtool(r, statistic="pvalue")
  fdr.out
}

##' Export differential correlations of comparison of two correlation matrices
##' 
##' @title Export differential correlations between two conditions
##' @param output.file can specify file name of the results exported
##' @param data1 data matrix under condition 1
##' @param data2 data matrix under condition 2
##' @param method c("pearson", "spearman", "kendall")
##' @param p.adjust.methods c("local", holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
##' @param threshold a threshold of significance levels of differential correlation
##' @return a text file
##' @examples
##' data(AraMetRoots)
##' AraMetRoots[AraMetRoots==0] <- NA
##' AraMetRootsImp <- completeObs(pca(log2(AraMetRoots), nPcs=3, method="ppca"))
##' comp.2.cc.fdr(output.file="res.txt", AraMetRootsImp[,1:17], method="spearman", 
##'               AraMetRootsImp[,18:37], threshold=0.05)
##' @references
##' Fukushima, A. Gene (2013) 518, 209-214
##' @author Atsushi Fukushima
##' @export
comp.2.cc.fdr <- function (output.file = "res.txt", data1, data2, method = "pearson", 
                           p.adjust.methods = "local", threshold = 0.05) 
{
  cc1 <- cor(t(data1), method = method)
  cc2 <- cor(t(data2), method = method)
  ccc1 <- as.vector(cc1[lower.tri(cc1)])
  ccc2 <- as.vector(cc2[lower.tri(cc2)])
  n1 <- ncol(data1)
  n2 <- ncol(data2)
  n <- nrow(data1)
  N <- n * (n - 1)/2
  p1 <- rep(1, N)
  p2 <- rep(1, N)
  pdiff <- rep(1, N)
  diff <- rep(1, N)
  mol.names <- rownames(cc1)
  p1 <- cor2.test(n1, ccc1)
  p2 <- cor2.test(n2, ccc2)
  pdiff <- compcorr(n1, ccc1, n2, ccc2)$pval
  diff <- ccc1 - ccc2
  pdiff[(is.na(pdiff)) == TRUE] <- 1
  if (p.adjust.methods == "local") {
    p1.lfdr <- get.lfdr(p1)$lfdr
    p2.lfdr <- get.lfdr(p2)$lfdr
    pdiff.lfdr <- get.lfdr(pdiff)$lfdr
  } else if (p.adjust.methods == "BH" | p.adjust.methods == "bh") {
    p1.lfdr <- p.adjust(p1, method=p.adjust.methods)
    p2.lfdr <- p.adjust(p2, method=p.adjust.methods)
    pdiff.lfdr <- p.adjust(pdiff, method=p.adjust.methods)
  } else {
    p1.lfdr <- rep("not adjusted", N)
    p2.lfdr <- rep("not adjusted", N)
    pdiff.lfdr <- rep("not adjusted", N)
  }
  
  ## generates combination
  myindex <- which((lower.tri(cc1))==TRUE, arr.ind=TRUE)
  mol.names1 <- mol.names[myindex[,2]]
  mol.names2 <- mol.names[myindex[,1]]
  
  ## threshold
  fin.ind <- pdiff.lfdr<threshold
  ## concat
  res <- cbind(mol.names1[fin.ind], mol.names2[fin.ind], ccc1[fin.ind], p1[fin.ind],
               ccc2[fin.ind], p2[fin.ind],
               pdiff[fin.ind], diff[fin.ind], p1.lfdr[fin.ind], p2.lfdr[fin.ind], pdiff.lfdr[fin.ind])
  head <- c("molecule X", "molecule Y", "r1", "p1",
            "r2", "p2", "p (difference)", "(r1-r2)",
            "lfdr (in cond. 1)", "lfdr (in cond. 2)",
            "lfdr (difference)")
  res <- rbind(head, res)
  write.table(res, file = output.file, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
}




##' Additional distance functions Correlation distance (1-r)
##'
##'
##' These functions were originally from 'hybridHclust' package.
##' We modified the functions slightly. See also the reference manual in detail.
##'
##' @title Additional distance functions correlation distance (1-r)
##' @param data  a data matrix ([data.frame object] row: metabolites,
##' col: samples or replicates)
##' @param methods a character string indicating which correlation coefficient is to be calculated. One of "pearson" (default), "spearman", or "kendall" can be abbreviated.
##' @param absolute TRUE means that absolute value of the correlation coefficient is used (Default: FALSE).
##' @return the resulting correlation matrix
##' @examples
##' cor.dist(as.matrix(t(iris[,1:4])))
##' @author Atsushi Fukushima
##' @export
cor.dist <- function(data, methods="pearson", absolute=FALSE) {
	my.dist <- cor(t(data), method=methods)
	if (absolute) {
          return(1-abs(my.dist))
        } else {
          return(1-my.dist)
        }
}

##' Additional distance functions correlation distance (uncentered)
##'
##'
##' These functions were originally from 'hybridHclust' package.
##' We modified the functions slightly. See also the reference manual in detail.
##'
##' @title Additional distance functions correlation distance (uncentered) 
##' @param data  a data matrix ([data.frame object] row: metabolites,
##' col: samples or replicates)
##' @param i i-th row of data
##' @param absolute TRUE means that absolute value of the correlation coefficient is used (Default: FALSE).
##' @return the resulting correlation matrix
##' @examples
##' uncent.cor2dist(as.matrix(t(iris[,1:4])), 1) ## NOT RUN!
##' @author Atsushi Fukushima
##' @export
uncent.cor2dist <- function(data, i, absolute=FALSE) {
# Calculates uncentered correlation distance
# between row i of data and all other rows of data.
  myfun <- function(y) {
    data1 <- data[i,]
    data1[is.na(y)] <- NA
    y[is.na(data1)] <- NA
    sum(y * data1 / sqrt( sum(y * y, na.rm=TRUE) *
                         sum(data1 * data1, na.rm=TRUE) ), na.rm=TRUE)
  }
  if(absolute) {
    return(1-abs(apply(data, 1, myfun)))
  }
  else {
    return(1-apply(data, 1, myfun))
  }
}
##' Calculating all pairwise distances using correlation distance
##'
##'
##' These functions were originally from 'hybridHclust' package.
##' We modified the functions slightly. See also the reference manual in detail.
##'
##' @title Calculating all pairwise distances using correlation distance
##' @param data  a data matrix ([data.frame object] row: metabolites,
##' col: samples or replicates)
##' @param absolute TRUE means that absolute value of the correlation coefficient is used (Default: FALSE).
##' @return the resulting correlation matrix
##' @examples
##' uncent.cordist(as.matrix(t(iris[,1:4]))) ## NOT RUN!
##' @author Atsushi Fukushima
##' @export
uncent.cordist <- function(data, absolute=FALSE) {
  m <- nrow(data)
  n <- ncol(data)
  my.dist <- matrix(0, m, n)
  for (i in 1:m) {
    my.dist[, i] <- uncent.cor2dist(data, i, absolute)
  }
  return(my.dist)  # value of "absolute" already used above
}


##' Cluster molecules
##'
##' @title Hierarchical clustering of molecules
##' @param data matrix or data frame
##' @param method c("pearson", "spearman", "kendall", "euclidean",
##' "maximum", "manhattan", "canberra", "binary", or "minkowski")
##' @param linkage c("average", "ward", "single", "complete", "mcquitty",
##' "median", "centroid")
##' @param absolute if TRUE, then 1-|COR| else 1-COR, default is FALSE
##' @return an object of class 'hclust'
##' @examples
##' cluster.molecule(as.matrix(t(iris[,1:4])), "pearson", "average")
##' @author Atsushi Fukushima
##' @export
cluster.molecule <- function(data, method="pearson", linkage="average", absolute=FALSE) {
  ## check input form
  if ( !is.matrix(data) & !is.data.frame(data) ) {
    stop("data must be a matrix or data frame.")
  }

  if ( match(method, c("pearson", "spearman", "kendall")) ) {
    dist <- as.dist( cor.dist(data, methods=method, absolute=absolute) )
  } else if ( match(method, c("euclidean","maximum","manhattan",
                               "canberra","binary","minkowski")) ) {
    dist <- dist(data, method=method)
  } else if ( match(method, c("uncentered")) ) {
    dist <- as.dist( uncent.cordist(data) )
  } else {
    stop("method must be a similarity or distance measure as described in cor(), dist(), or uncent.cordist.")
  }
  hclust(d=dist, method=linkage)
}




##' Get eigen molecule
##'
##' @title Get eigen molecule
##' @param data  a data matrix ([data.frame object] row: molecules,
##' col: samples or replicates)
##' @param groups a vector of group memberships as returned by cutree
##' @param whichgroups a vector of group numbers to examine
##' @param methods c("svd", "nipals", "rnipals", "bpca", "ppca"). See also pca() function in pcaMethods package
##' @param n top n principal components
##' @return the resulting vector.
##' @examples
##' library(pcaMethods)
##' data(golub, package = "multtest")
##' hc.mol1 <- cluster.molecule(golub[1:100, 1:27], "pearson", "average")
##' g1 <- cutree(hc.mol1, h=0.6)
##' res1 <- get.eigen.molecule(golub[1:100,], g1)
##' @author Atsushi Fukushima
##' @export
get.eigen.molecule <- function(data, groups, whichgroups=NULL, methods="svd", n=10){
  ## grouping
  gs <- data.frame(group=unique(groups), size=as.vector(by(data, groups, nrow)))
  # which groups are we looking at?
  if (is.null(whichgroups)) {
     # by default we look at all groups where n > 10
     whichgroups <- gs$group[gs$size >= n]
  }
  # set up vectors to hold the return information
  group <- vector()
  N <- vector()
  eigen.molecules <- list()
  mean.corr <- vector()
  # index
  no.res <- 1
  for (i in whichgroups) {
    print(no.res)
    ## calculate size of the group
    d.length <- gs$size[gs$group == i]
    ## dividing an orginal matrix into submatrix
    submat <- data[groups==i,]
    ## correlation matrix
    c1 <- cor(t(submat))
    ## PCA
    mydata <- submat
    pc <- pca(t(mydata), method=methods, nPcs=n)
    eigen.molecules[[no.res]] <- as.vector(loadings(pc)[1:n,1])    
    group[no.res] <- i
    N[no.res] <- nrow(data[groups==i,])
    mean.corr[no.res] <- mean(c1[lower.tri(c1)]) 
    no.res <- no.res + 1
  }
  df <- list("group"=group, "N"=N, "mean.corr"=mean.corr, "eigen.molecules"=eigen.molecules)
}



##' Generating graph from data matrix
##'
##' @title Generating graph from data matrix
##' @param data data matrix or data frame
##' @param method c("Pearson", "Spearman", "Kendall")
##' @param cor.thr a threshold of correlation coefficient (default: r >= 0.6)
##' @param neg.flag flag where uses or not negative correlations
##' @param node.col specifies color of nodes in a graph (default: red)
##' @param node.size specifies size of nodes in a graph (default: 7)
##' @param edge.col specifies color of edges in a graph (default: blue)
##' @param edge.width specifies width of edges in a graph (default: 3)
##' @return igraph object
##' @examples
##' library(igraph)
##' mat <- matrix(runif(100), nr=10)
##' rownames(mat) <- as.character(1:10)
##' generate_g(mat)
##' @author Atsushi Fukushima
##' @export
generate_g <- function(data, method="pearson", cor.thr=0.6, neg.flag=1,
                       node.col="red", node.size=7, edge.col="blue", edge.width=3) {
  # calc correlation matrix
  data.cor <- cor(t(data), method=method)
  # association network
  if (neg.flag==1) {
    g <- simplify( graph.adjacency( (data.cor >= cor.thr |
                                     data.cor < -cor.thr), mode="undirected") )
  } else {
    g <- simplify( graph.adjacency( (data.cor >= cor.thr), mode="undirected") )
  }
  V(g)$size <- node.size
  V(g)$name <- row.names(data)
  V(g)$color <- node.col
  V(g)$label <- row.names(data)
  E(g)$color <- edge.col
  E(g)$width <- edge.width
  g
}


##' Getting graph from eigengene module list
##'
##' @title Getting graph from eigengene module list
##' @param eigen.list the resulting vector from get.eigen.molecule
##' @param label a label of module extracted (default: "Module")
##' @return igraph object
##' @examples
##' library(pcaMethods)
##' library(igraph)
##' data(golub, package = "multtest")
##' hc.mol1 <- cluster.molecule(golub[, 1:27], "pearson", "average")
##' g1 <- cutree(hc.mol1, h=0.4)
##' res1 <- get.eigen.molecule(golub, g1)
##' g1.eigen <- get.eigen.molecule.graph(res1)
##' @author Atsushi Fukushima
##' @export
get.eigen.molecule.graph <- function(eigen.list, label="Module") {
  res <- eigen.list
  mat <- t(matrix(unlist(res$eigen.molecules),
                       nrow=length(res$eigen.molecules[[1]])))
  rownames(mat) <- paste(label, 1:length(res$eigen.molecules), sep="")

  g <- generate_g(mat)
}


##' Writing modules into a text file
##'
##' @title Writing modules into a text file
##' @param cutree.res the result of cutree function
##' @param mod.list the result of get.eigen.molecule
##' @param outfile file name of output
##' @return a text file
##' @examples
##' data(golub, package = "multtest")
##' hc.mol1 <- cluster.molecule(golub[, 1:27], "pearson", "average")
##' g1 <- cutree(hc.mol1, h=0.4)
##' res1 <- get.eigen.molecule(golub, g1)
##' write.modules(g1, res1)
##' @author Atsushi Fukushima
##' @export
write.modules <- function(cutree.res,  mod.list, outfile="module_list.txt") {
  write("Molecule number\tFeature name (probeset or metabolite)",
        file=outfile, append=FALSE)
  for (i in mod.list$group) {
    targets <- cutree.res[cutree.res == i]
    res <- cbind(i, names(targets))
    write.table(res, file=outfile, quote=FALSE, row.names=FALSE,
                col.names=FALSE, sep="\t", append=TRUE)
  }
}

##' Get minimum and maximum
##'
##' @title Get minimum and maximum
##' @param d data matrix or data frame
##' @return list object of minimum value or maximum value in a data
##' @examples
##' get.min.max(iris[,1:2])
##' @author Atsushi Fukushima
##' @export
get.min.max <- function(d) {
  maxloc  <- which.max(as.vector(as.matrix(d)))/nrow(d)
  max.col <- ceiling(maxloc)
  max.row <- round((maxloc-trunc(maxloc)) * nrow(d))
  if (maxloc%%1 == 0) {
    max.row <- nrow(d)
  }
  
  minloc  <- which.min(as.vector(as.matrix(d)))/nrow(d)
  min.col <- ceiling(minloc)
  min.row <- round((minloc-trunc(minloc)) * nrow(d))
  if (minloc%%1 == 0) {
    min.row <- nrow(d)
  }
  
  list(max = d[max.row,max.col], min = d[min.row,min.col], max.i = max.row)
}



##' Plot cluster molecules
##'
##' @title Plot cluster molecules
##' @param data data matrix or data frame
##' @param groups a vector of group memberships as returned by cutree
##' @param group.no the group number to be plotted
##' @param title a title for the graph
##' @param ylim a vector indicating the upper and lower limit for the y-axis
##' @param order whether or not to order the columns of the data matrix
##' @param scale.center unless NULL, each row is scaled using scale
##' @param scale.scale unless NULL, each row is scaled using scale.
##' @param frame the color of the frame that is drawn as the background of the plot
##' @param col If NULL, all genes will be drawn in the default color (blue).  
##'  If the text "random" is given, then a set of colors will be generated by 
##' @param bottom.mar The size of the bottom margin of the plots as sent in par(mar=c(...))
##' @param xlab a lalel of x axis (defalt: "Samples")
##' @param ylab a lalel of y axis (defalt: "Relative abundance")
##' @return a graph
##' @references
##' this function was originally from Watson M (2005) BMC Bioinformatics 7:509
##' @examples
##' 
##'
##' library(pcaMethods)
##' data(golub, package = "multtest")
##' hc.mol1 <- cluster.molecule(golub[, 1:27], "pearson", "average")
##' g1 <- cutree(hc.mol1, h=0.4)
##' plotClusterMolecules(golub[,1:27], g1, 3)
##' @author Atsushi Fukushima
##' @export
plotClusterMolecules <-
function(data, groups=NULL, group.no=NULL, title=NULL, ylim=NULL, order=NULL,
         scale.center=FALSE, scale.scale=FALSE, frame="white", col=NULL,
         bottom.mar=5, xlab="Samples", ylab="Relative abundance") {

	## data is a data frame
	## groups is a vector as produced by cutree
	## group.no is the group you want to plot

	if (!is.null(groups) & !is.null(group.no)) {
		data <- as.matrix(data[groups==group.no,])
	}
	else {
		data <- as.matrix(data)
	}

	if (!is.null(scale.center) | !is.null(scale.scale)) {
		# scale the rows of the matrix
		data <- t(scale(t(data), center=scale.center, scale=scale.scale))
	}

	red <- NULL

	if(!is.null(order)) {

		# convert to data.frame
		data <- as.data.frame(data)

		if (order == "average") {
			d1 <- order(mean(data))
			red <- mean(data)
			red <- t(red)[,d1]
		}
		else {
			d1 <- order(t(data[rownames(data) == order,]))
			red <- data[rownames(data) == order, d1]
		}
		data <- data[, d1]
	}

	## find the min and max number
	mm <- get.min.max(data)

	if (is.null(ylim)) {
		ylim <- range(mm$min,mm$max)
	}
	
	if(is.null(col)) {
		# if col is null, we want all blue
		col <- rep("black", nrow(data))
	} else if (col == "random") {
		# we need to generate colors
		col<-rainbow(nrow(data))
	} else if (length(col) == 1) {
    # if given a single color, we want all
    col <- rep(col, nrow(data))
  } else {
 	  # we need to generate colors
		col<-rainbow(nrow(data))
  }

	## plot the largest gene
	par(ps=8, las=2, mar=c(bottom.mar,4,3,2)) # lwd=2
	plot(xy.coords(1:ncol(data), data[mm$max.i,]), type="l", col=col[mm$max.i], ylim=ylim, xlab=xlab, ylab=ylab, xaxt="n")
	axis(1, at=1:ncol(data), colnames(data))

	## plot the grey background
	rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = frame)

	## plot the rest
	for (i in 1:nrow(data)) {
		 lines(xy.coords(1:ncol(data),data[i,]), col=col[i])
	}

	## plot the red gene
	if(!is.null(red)) {
		lines(xy.coords(1:ncol(data),red), col="red", lty="dashed")
	}
	
	if(!is.null(title)) {
	   title(title)
	}
}



##' Plot DiffCorr group
##'
##' @title Plot DiffCorr group
##' @param data a data matrix or data frame
##' @param groups1 a vector of row group membership as produced by cutree under condition 1
##' @param groups2 a vector of row group membership as produced by cutree under condition 2
##' @param group1.no the group number to be plotted (condition 1)
##' @param group2.no the group number to be plotted (condition 2)
##' @param g1 a vector describing the columns of the data belonging to condition 1
##' @param g2 a vector describing the columns of the data belonging to condition 2
##' @param g1.order whether or not to order the columns of the data
##' matrix for condition 1.  If "average", then the columns are ordered
##' by the average expression value.  If the name of a gene (row),
##' then the columns are orderd according to the expression levels of
##' that gene.  If NULL, columns remain in their original order.  
##' @param g2.order See g1.order
##' @param title1 A title for the left hand graph
##' @param title2 A title for the right hand graph
##' @param ... other parameters to be passed to this function
##' @return a graph
##' @examples
##' library(pcaMethods)
##' data(golub, package = "multtest")
##' hc.mol1 <- cluster.molecule(golub[, 1:27], "pearson", "average")
##' hc.mol2 <- cluster.molecule(golub[, 28:38], "pearson", "average")
##' g1 <- cutree(hc.mol1, h=0.4)
##' g2 <- cutree(hc.mol2, h=0.4)
##' ##
##' plotDiffCorrGroup(golub, g1, g2, 21, 24, 1:27, 28:38,
##'                    scale.center=TRUE, scale.scale=TRUE,
##'                    ylim=c(-5,5))
##' @author Atsushi Fukushima
##' @export
plotDiffCorrGroup <-
function(data, groups1=NULL, groups2=NULL, group1.no=NULL, group2.no=NULL, g1, g2, g1.order=NULL, g2.order=NULL, title1=NULL, title2=NULL,...) {
	if ((!is.null(groups1) & !is.null(groups2)) & (!is.null(group1.no) & (!is.null(group2.no))) ) {
		data1 <- as.matrix(data[groups1==group1.no,])
		data2 <- as.matrix(data[groups2==group2.no,])                
	}
	else {
		data <- as.matrix(data)
	}
        ##
	g1.data  <- data1[, g1]
	g2.data  <- data2[, g2]
        ##
	split.screen(c(1,2))
	screen(1)
	plotClusterMolecules(data = g1.data, order = g1.order, title = title1,...)
	screen(2)
	plotClusterMolecules(data = g2.data, order = g2.order, title = title2,...)
	close.screen(all.screens = TRUE)
}




