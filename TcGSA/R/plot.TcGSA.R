#'Plot a Gene Set Trends Heatmap.
#'
#'This function plots a gene sets dynamic trends heatmap.
#'
#'On the heatmap, each line corresponds to a gene set, and each column to a
#'timepoint.
#'
#'If \code{expr} is a matrix or a dataframe, then the "original" data are
#'plotted.  On the other hand, if \code{expr} is a list returned in the
#'\code{'Estimations'} element of \code{\link{TcGSA.LR}}, then it is those
#'"estimations" made by the \code{\link{TcGSA.LR}} function that are plotted.
#'
#'If \code{descript} is \code{FALSE}, the second element of \code{margins} can
#'be reduced (for instance use \code{margins = c(5, 10)}), as there is not so
#'much need for space in order to display only the gene set names, without
#'their description.
#'
#'If there is a large number of significant gene sets, the hierarchical clustering
#'step repeated for each of them can take a few minutes. To speed things up 
#'(especially) when playing with the ploting parameters for having a nice plot,
#'one can run the \code{clustTrend} function beforehand, and plug its results 
#'in the \code{plot.TcGSA} function via the \code{clust_trends} argument.
#'
#'
#'@method plot TcGSA
#'
#'@param x 
#'an object of class'\code{TcGSA}'.
#'
#'@param threshold 
#'the threshold at which the FDR or the FWER should be
#'controlled.
#'
#'@param myproc 
#'a vector of character strings containing the names of the
#'multiple testing procedures for which adjusted p-values are to be computed.
#'This vector should include any of the following: "\code{Bonferroni}",
#'"\code{Holm}", "\code{Hochberg}", "\code{SidakSS}", "\code{SidakSD}",
#'"\code{BH}", "\code{BY}", "\code{ABH}", "\code{TSBH}" or "\code{none}".  
#'"\code{none}" indicates no adjustement for multiple testing. See
#'\code{\link[multtest:mt.rawp2adjp]{mt.rawp2adjp}} for details.  Default is
#'"\code{BY}", the Benjamini & Yekutieli (2001) step-up FDR-controlling
#'procedure (general dependency structures).  In order to control the FWER(in
#'case of an analysis that is more a hypothesis confirmation than an
#'exploration of the expression data), we recommand to use "\code{Holm}", the
#'Holm (1979) step-down adjusted p-values for strong control of the FWER.
#'
#'@param nbsimu_pval 
#'the number of observations under the null distribution to
#'be generated in order to compute the p-values. Default is \code{1e+06}.
#'
#'@param expr 
#'either a matrix or dataframe of gene expression upon which
#'dynamics are to be calculated, or a list of gene sets estimation of gene
#'expression.  In the case of a matrix or dataframe, its dimension are \eqn{n}
#'x \eqn{p}, with the \eqn{p} sample in column and the \eqn{n} genes in row.
#'In the case of a list, its length should correspond to the number of gene
#'sets under scrutiny and each element should be an 3 dimension array of
#'estimated gene expression, such as for the list returned in the
#'\code{'Estimations'} element of \code{\link{TcGSA.LR}}.  See details.
#'
#'@param Subject_ID 
#'a factor of length \eqn{p} that is in the same order as the
#'columns of \code{expr} (when it is a dataframe) and that contains the patient
#'identifier of each sample. Ignored if \code{expr} is a list of estimations.
#TODO See Details.
#'
#'@param TimePoint 
#'a numeric vector or a factor of length \eqn{p} that is in
#'the same order as \code{Subject_ID} and the columns of \code{expr} (when it
#'is a dataframe), and that contains the time points at which gene expression
#'was measured. Ignored if \code{expr} is a list of estimations.
#TODO See Details.
#'
#'@param baseline 
#'a character string which is the value of \code{TimePoint}
#'used as baseline.  See Details.
#'
#'@param only.signif 
#'logical flag for plotting only the significant gene sets.
#'If \code{FALSE}, all the gene sets from the \bold{gmt} object contained in
#'\code{x} are plotted.  Default is \code{TRUE}.
#'
#'@param group.var
#'in the case of several treatment` groups, this is a factor of
#'length \eqn{p} that is in the same order as \code{Timepoint},
#'\code{Subject_ID}, \code{sample_name} and the columns of \code{expr}.  It
#'indicates to which treatment group each sample belongs to.  Default is
#'\code{NULL}, which means that there is only one treatment group.  See
#'Details.
#'
#'@param Group_ID_paired
#' a character vector of length \eqn{p} that is in the
#'same order as \code{Timepoint}, \code{Subject_ID}, \code{sample_name},
#'\code{group.var} and the columns of \code{expr}.  This argument must not be
#'\code{NULL} in the case of a paired analysis, and must be \code{NULL}
#'otherwise.  Default is \code{NULL}.  See Details.
#'
#'@param ref 
#'the group which is used as reference in the case of several
#'treatment groups.  Default is \code{NULL}, which means that reference is the
#'first group in alphabetical order of the labels of \code{group.var}.  See
#'Details.
#'
#'@param group_of_interest
#'the group of interest, for which dynamics are to be
#'computed in the case of several treatment groups.  Default is \code{NULL},
#'which means that group of interest is the second group in alphabetical order
#'of the labels of \code{group.var}.  See Details.%% ~~Describe
#'\code{group_of_interest} here~~
#'
#'@param ranking
#'a logical flag. If \code{TRUE}, the gene set trends are not hierarchicaly classified, but
#'ordered by decreasing Likelihood ratios. Default is \code{FALSE}.
#'
#'@param FUNcluster 
#'the clustering function used to agglomerate genes in
#'trends.  Default is \code{NULL}, in which a hierachical clustering is
#'performed via the function \code{\link[cluster:agnes]{agnes}}, using the
#'metric \code{clustering_metric} and the method \code{clustering_method}.  See
#'\code{\link[cluster:clusGap]{clusGap}}
#'
#'@param clustering_metric 
#'character string specifying the metric to be used
#'for calculating dissimilarities between observations in the hierarchical
#'clustering when \code{FUNcluster} is \code{NULL}.  The currently available
#'options are \code{"euclidean"} and \code{"manhattan"}.  Default is
#'\code{"euclidean"}.  See \code{\link[cluster:agnes]{agnes}}.  Also, a \code{"sts"} option 
#'is available in TcGSA.  It implements the 'Short Time Series' distance 
#'[Moller-Levet et al., Fuzzy CLustering of short time series and unevenly distributed 
#'sampling points, \emph{Advances in Intelligent Data Analysis V}:330-340 Springer, 2003]
#'designed specifically for clustering time series.
#'
#'@param clustering_method 
#'character string defining the agglomerative method
#'to be used in the hierarchical clustering when \code{FUNcluster} is
#'\code{NULL}.  The six methods implemented are \code{"average"} ([unweighted
#'pair-]group average method, UPGMA), \code{"single"} (single linkage),
#'\code{"complete"} (complete linkage), \code{"ward"} (Ward's method),
#'\code{"weighted"} (weighted average linkage).  Default is \code{"ward"}.  See
#'\code{\link[cluster:agnes]{agnes}}.
#'
#'@param B 
#'integer specifying the number of Monte Carlo ("bootstrap") samples
#'used to compute the gap statistics.  Default is \code{500}.  See
#'\code{\link[cluster:clusGap]{clusGap}}.
#'
#'@param max_trends 
#'integer specifying the maximum number of different clusters
#'to be tested.  Default is \code{4}.
#'
#'@param aggreg.fun 
#'a character string such as \code{"mean"}, \code{"median"}
#'or the name of any other statistics function defined that returns a single
#'numeric value.  It specifies the function used to aggregate the observations
#'before the clustering.  Default is to \code{median}.  Default is
#'\code{"median"}.
#'
#'@param methodOptiClust 
#'character string indicating how the "optimal"" number
#'of clusters is computed from the gap statistics and their standard
#'deviations. Possible values are \code{"globalmax"}, \code{"firstmax"},
#'\code{"Tibs2001SEmax"}, \code{"firstSEmax"} and \code{"globalSEmax"}.
#'Default is \code{"firstSEmax"}.  See \code{'method'} in
#'\code{\link[cluster:clusGap]{clusGap}}, Details and \emph{Tibshirani et al.,
#'2001} in References.
#'
#'@param indiv 
#'a character string indicating by which unit observations are
#'aggregated (through \code{aggreg.fun}) before the clustering.  Possible
#'values are \code{"genes"} or \code{"patients"}.  Default is \code{"genes"}.
#'
#'@param verbose 
#'logical flag enabling verbose messages to track the computing
#'status of the function.  Default is \code{TRUE}.
#'
#'@param clust_trends 
#'object of class \bold{\link{ClusteredTrends}} containing
#'already computed trends for the plotted gene sets.  Default is \code{NULL}.
#'
#'@param N_clusters 
#'an integer that is the number of clusters in which the
#'dynamics should be regrouped. The cutoff of the clustering tree is
#'automatically calculated accordingly.  Default is \code{NULL}, in which case
#'the dendrogram is not cut and no clusters are identified.
#'
#'@param myclusters 
#'a character vector of colors for predefined clusters of the
#'represented genesets, with as many levels as the value of \code{N_clusters}.
#'Default is \code{NULL}, in which case the clusters are automatically
#'identified and colored via the \code{\link[stats:cutree]{cutree}} function and the
#'\code{N_clusters} argument only.
#'
#'@param label.clusters 
#'if \code{N_clusters} is not \code{NULL}, a character
#'vector of length \code{N_clusterss}.  Default is \code{NULL}, in which case
#'if \code{N_clusters} is not \code{NULL}, clusters are simply labelled with
#'numbers.
#'
#'@param prev_rowCL 
#'a \bold{hclust} object, such as the one return by the
#'present plotting funstion (see Value) for instance.  If not \code{NULL}, no
#'clustering is calculated by the present plotting function and this tree is
#'used to represent the gene sets dynamics.  Default is \code{NULL}.
#'
#'@param descript 
#'logical flag indicating that the description of the gene sets
#'should appear after their name on the right side of the plot if \code{TRUE}.
#'Default is \code{TRUE}.  See Details.
#'
#'@param plot 
#'logical flag indicating wether the heatmap should be plotted or
#'not.  Default is \code{TRUE}.
#'
#'@param color.vec 
#'a character strings vector used to define the color
#'\link[grDevices:palette]{palette} used in the plot.  Default is
#'\code{c("#D73027", "#FC8D59","lightyellow", "#91BFDB", "#4575B4")}.
#'
#'@param legend.breaks 
#'a numeric vector indicating the splitting points for
#'coloring.  Default is \code{NULL}, in which case the break points will be
#'spaced equally and symetrically about 0.
#'
#'@param label.column 
#'a vector of character strings with the labels to be
#'displayed for the columns (i.e. the time points).  Default is \code{NULL}.
#'
#'@param time_unit 
#'the time unit to be displayed (such as \code{"Y"},
#'\code{"M"}, \code{"W"}, \code{"D"}, \code{"H"}, etc) next to the values of
#'\code{TimePoint} in the columns labels when \code{label.column} is
#'\code{NULL}.  Default is \code{""}.
#'
#'@param cex.label.row 
#'a numerical value giving the amount by which row labels
#'text should be magnified relative to the default \code{1}.
#'
#'@param cex.label.column 
#'a numerical value giving the amount by which column
#'labels text should be magnified relative to the default \code{1}.
#'
#'@param margins 
#'numeric vector of length 2 containing the margins (see
#'\code{\link{par}(mar= *))} for column and row names, respectively.  Default
#'is \code{c(15, 100)}.  See Details.
#'
#'@param heatKey.size 
#'the size of the color key for the heatmap fill.  Default
#'is \code{1}.
#'
#'@param dendrogram.size 
#'the horizontal size of the dendrogram.  Default is \code{1}
#'
#'@param heatmap.height 
#'the height of the heatmap.  Default is \code{1}
#'
#'@param heatmap.width 
#'the width of the heatmap.  Default is \code{1}
#'
#'@param cex.clusterKey 
#'a numerical value giving the amount by which the
#'clusters legend text should be magnified relative to the default \code{1},
#'when \code{N_clusters} is not \code{NULL}.
#'
#'@param cex.main 
#'a numerical value giving the amount by which title text
#'should be magnified relative to the default \code{1}.
#'
#'@param horiz.clusterKey 
#'a logical flag; if \code{TRUE}, set the legend for
#'clusters horizontally rather than vertically.  Only used if the
#'\code{N_clusters} argument is not \code{NULL}.  Default is \code{TRUE}.
#'@param main 
#'a character string for an optionnal title.  Default is \code{NULL}.
#'
#'@param subtitle
#' a character string for an optionnal subtitle.  Default is \code{NULL}.
#'
#'@param \dots 
#'other parameters to be passed through to plotting functions.
#'
#'@return An object of class \bold{\link[stats:hclust]{hclust}} which describes the tree
#'produced by the clustering process.  The object is a list with components:
#'\itemize{
#'\item \code{merge} an \eqn{n-1} by \eqn{2} matrix.  Row \eqn{i} of
#'\code{merge} describes the merging of clusters at step i of the clustering.
#'If an element \eqn{j} in the row is negative, then observation -\eqn{j} was
#'merged at this stage.  If \eqn{j} is positive then the merge was with the
#'cluster formed at the (earlier) stage \eqn{j} of the algorithm.  Thus
#'negative entries in merge indicate agglomerations of singletons, and positive
#'entries indicate agglomerations of non-singletons.
#'\item \code{height} a set of \eqn{n-1} real values (non-decreasing for
#'ultrametric trees).  The clustering height: that is, the value of the
#'criterion associated with the Ward clustering method.
#'\item \code{order} a vector giving the permutation of the original
#'observations suitable for plotting, in the sense that a cluster plot using
#'this ordering and matrix merge will not have crossings of the branches.
#'\item \code{labels} the gene set trends name.
#'\item \code{call} the call which produced the result clustering:
#'\cr\code{hclust(d = dist(map2heat, method = "euclidean"), method = "ward.D2")}
#'\item \code{method} "ward.D2", as it is the clustering method that has been used
#'for clustering the gene set trends.
#'\item \code{dist.method} "euclidean", as it is the distance that has been used
#'for clustering the gene set trends.
#'\item \code{legend.breaks} a numeric vector giving the splitting points used
#'for coloring the heatmap.  If \code{plot} is \code{FALSE}, then it is
#'\code{NULL}.
#'\item \code{myclusters} a character vector of colors for the dynamic clusters
#'of the represented gene set trends, with as many levels as the value of
#'\code{N_clusters}.  If no dynamic clusters were represented, than this is
#'\code{NULL}.
#'\item \code{ddr} a \bold{dendrogram} object with the reordering used for the
#'heatmap.  See \code{\link[gplots:heatmap.2]{heatmap.2}}.
#'\item geneset.names character vector with the names of the gene sets
#'used in the heatmap.
#'\item \code{clust.trends} a \bold{\link{ClusteredTrends}} object.
#'\item \code{clustersExport} a data frame with 2 variables containing the two
#'following variables : \itemize{ \item \code{GeneSet}: the gene set trends
#'clustered.  \item \code{Cluster}: the dynamic cluster they belong to.  } The
#'data frame is order by the variable \code{Cluster}.
#'\item \code{data_plotted}: the data matrix represented by the heatmap
#'}
#'
#'@author Boris P. Hejblum
#'
#'@seealso \code{\link[gplots:heatmap.2]{heatmap.2}}, \code{\link{TcGSA.LR}},
#'\code{\link[stats:hclust]{hclust}}
#'
#'@references Hejblum BP, Skinner J, Thiebaut R, (2015) 
#'Time-Course Gene Set Analysis for Longitudinal Gene Expression Data. 
#'\emph{PLoS Computat Biol} 11(6): e1004310.
#'doi: 10.1371/journal.pcbi.1004310
#'
#'@import ggplot2
#'
#'@import reshape2
#'
#'@importFrom graphics legend
#'
#'@importFrom grDevices colorRampPalette hsv rainbow
#'
#'@importFrom stats as.dendrogram cutree dist hclust order.dendrogram quantile reorder
#'
#'@export
#'
#'@examples
#'
#'data(data_simu_TcGSA)
#'
#'tcgsa_sim_1grp <- TcGSA.LR(expr=expr_1grp, gmt=gmt_sim, design=design, 
#'                           subject_name="Patient_ID", time_name="TimePoint",
#'                           time_func="linear", crossedRandom=FALSE)
#'summary(tcgsa_sim_1grp)
#'
#'plot(x=tcgsa_sim_1grp, expr=tcgsa_sim_1grp$Estimations, 
#'     Subject_ID=design$Patient_ID, TimePoint=design$TimePoint,
#'     baseline=1, 
#'     B=100,
#'     time_unit="H",
#'     dendrogram.size=0.4, heatmap.width=0.8, heatmap.height=2, cex.main=0.7
#'     )
#'
#'\dontrun{                
#'tcgsa_sim_2grp <- TcGSA.LR(expr=expr_2grp, gmt=gmt_sim, design=design, 
#'                           subject_name="Patient_ID", time_name="TimePoint",
#'                           time_func="linear", crossedRandom=FALSE, 
#'                           group_name="group.var")
#'summary(tcgsa_sim_2grp)                             
#'plot(x=tcgsa_sim_2grp, expr=expr_2grp, 
#'     Subject_ID=design$Patient_ID, TimePoint=design$TimePoint,
#'     B=100,
#'     time_unit="H",
#'     )
#'}
#'
#'
plot.TcGSA <-
	function(x, threshold=0.05, myproc="BY", nbsimu_pval=1e+06, 
			 expr, Subject_ID, TimePoint, 
			 baseline=NULL, only.signif=TRUE,
			 group.var=NULL, Group_ID_paired=NULL, ref=NULL, group_of_interest=NULL,
			 ranking=FALSE,
			 FUNcluster=NULL, clustering_metric="euclidian", clustering_method="ward", B=500,
			 max_trends=4, aggreg.fun="median",
			 methodOptiClust = "firstSEmax",
			 indiv="genes",
			 verbose=TRUE,
			 clust_trends=NULL,
			 N_clusters=NULL, myclusters=NULL, label.clusters=NULL, prev_rowCL=NULL,
			 descript=TRUE, plot=TRUE,
			 color.vec=c("darkred", "#D73027", "#FC8D59", "snow", "#91BFDB", "#4575B4", "darkblue"),
			 legend.breaks=NULL,
			 label.column=NULL, time_unit="", 
			 cex.label.row=1, cex.label.column=1, margins=c(5, 25), heatKey.size=1, dendrogram.size=1, heatmap.height=1, heatmap.width=1,
			 cex.clusterKey=1, cex.main=1,
			 horiz.clusterKey=TRUE,
			 main=NULL, subtitle=NULL, 
			 ...){
		
		gmt <- x[["GeneSets_gmt"]]
		
		if(!is.null(baseline)){
			if(is.data.frame(expr) | is.matrix(expr)){
				if(!(baseline %in% unique(TimePoint))){
					stop("The 'baseline' value used is not one of the time points in 'TimePoint'...\n\n")
				}
			}
			else if(is.list(expr)){
				if(!(baseline %in% dimnames(expr[[1]])[[3]])){
					stop("The 'baseline' value used is not one of the time points in 'TimePoint'...\n\n")
				}
			}
		}
		
		if(is.null(main)){
			if(!is.null(subtitle)){
				mymain <-paste("Median trends", "\n", subtitle, sep="")
			}else{
				mymain <-paste("Median trends\n", "over all patients", sep="")
			}
		}else {
			if(!is.null(subtitle)){
				mymain <-paste(main, "\n", subtitle, sep="")
			}else{
				mymain <- main
			}
		}
		
		if(is.null(prev_rowCL)){
			clRows=TRUE
			signif <- multtest.TcGSA(x, threshold, myproc, nbsimu_pval)
			select <- which(signif$adj_pval<0.05)
			subtitle <- paste(subtitle, "\n", length(which(signif$adj_pval<0.05)), "/", length(signif$adj_pval), " gene sets significant with ", x[["func_form"]], " shape", sep="")
		}else{
			if(!is.null(prev_rowCL$ddr)){
				clRows=prev_rowCL$ddr
			}else{
				clRows=stats::as.dendrogram(prev_rowCL)
			}
			
			select <- match(prev_rowCL$geneset.names, gmt$geneset.names)
			if(length(which(is.na(select)))>0){
				select <- match(prev_rowCL$labels, gsub(": Undetermined","", paste(gmt$geneset.names, ": ", gmt$geneset.description, sep="")))
			}
			if(length(which(is.na(select)))>0){
				stop("Geneset names used in the previous clustering don't match the 'geneset.names' from the 'gmt' element of the 'x' argument")
			}
			
			if(is.null(myclusters) && !is.null(prev_rowCL$myclusters)){
				myclusters <- prev_rowCL$myclusters
			}  
		}
		
		if(only.signif & is.null(clustTrend)){
			if(!length(select)>0){
				stop("No gene sets significant")
			}
			gmt <- list("genesets"=gmt$genesets[select], "geneset.names"=gmt$geneset.names[select], "geneset.descriptions"=gmt$geneset.descriptions[select])
			class(gmt) = "GSA.genesets"
		}
		if(is.null(clust_trends)){
			if (length(which(!is.na(x$fit$LR)))<1){
				stop ("SERIOUS PROBLEM\n Was not able to compute any likelihood ratios...")
			}
			clust_trends <- clustTrend(tcgs=x, expr=expr, Subject_ID=Subject_ID, TimePoint=TimePoint, baseline=baseline, only.signif=only.signif,
									   group.var=group.var, Group_ID_paired=Group_ID_paired, ref=ref, group_of_interest=group_of_interest,
									   FUNcluster=FUNcluster, clustering_metric=clustering_metric, clustering_method=clustering_method, B=B,
									   max_trends=max_trends, aggreg.fun=aggreg.fun,
									   methodOptiClust = methodOptiClust,
									   indiv=indiv,
									   verbose=verbose
			)
		}else if(class(clust_trends)!="ClusteredTrends"){
			stop("The 'clust_trends' argument is not of the class 'ClusteredTrends', see the clustTrend function")
		}
		
		medoids2clust <- reshape2::acast(reshape2::melt(clust_trends[["ClustMeds"]], variable.name="Cluster", id.vars="TimePoint"),
										 formula="L1 + Cluster~ TimePoint", value.var="value")
		gsNames <- gsub("_.*$", "", rownames(medoids2clust))
		ncl <- gsub("^.*?_", "", rownames(medoids2clust))
		medoids2clust <- medoids2clust[,order(as.numeric(colnames(medoids2clust)))]
		
		if(ranking){
			trendsIndex <- match(gsNames, gmt$geneset.names)
			LRtrends <- x$fit$LR[trendsIndex]
			rank <- order(LRtrends, decreasing=TRUE)
			gsNames <- gsNames[rank]
			ncl <- ncl[rank]
			medoids2clust <- medoids2clust[rank, ]
			
			LR2quant <- x$fit$LR
			LR2quant[which(is.na(LR2quant))] <- 0
			percentiles <- stats::quantile(LR2quant, probs=seq(0.01, 1, 0.01))
			percTrends <- findInterval(LRtrends[rank], vec=percentiles)  
		}
		
		if(!descript){
			rownames(medoids2clust) <- paste(gsub("_", " ", rownames(medoids2clust)), clust_trends[[1]][gsNames], sep="/")
		}else{
			if(ranking){
				rownames(medoids2clust) <- paste(gsub(": Undetermined","", paste(gmt$geneset.names[match(gsNames, gmt$geneset.names)], 
																				 "[", percTrends, "th pctile]: ", 
																				 gmt$geneset.description[match(gsNames, gmt$geneset.names)], sep="")),
												 " ", ncl, "/", clust_trends[[1]][gsNames], sep="")
			}else{
				rownames(medoids2clust) <- paste(gsub(": Undetermined","", paste(gmt$geneset.names[match(gsNames, gmt$geneset.names)], ": ", 
																				 gmt$geneset.description[match(gsNames, gmt$geneset.names)], sep="")),
												 " ", ncl, "/", clust_trends[[1]][gsNames], sep="")
			}
		}
		
		map2heat <- medoids2clust
		# browser()
		# map2heat <- t(apply(X=medoids2clust, MARGIN=1, FUN=scale, center=FALSE)) # scales each rows
		# dimnames(map2heat) <- dimnames(medoids2clust)
		# map2heat <-  map2heat[grep("]:", rownames(map2heat)), ] # selects only annotated gene sets
		
		if(is.null(prev_rowCL)){
			hc <- stats::hclust(d=stats::dist(map2heat, method = "euclidean"), method="ward.D2")
			row_wt <- rowMeans(x=map2heat, na.rm = TRUE)
			ddr <- stats::reorder(x=stats::as.dendrogram(hc), wts=row_wt)
		}else{
			ddr <- prev_rowCL$ddr
			hc <- prev_rowCL
		}
		
		
		
		if(!is.null(N_clusters) && is.null(myclusters)){
			myclusters <- as.factor(stats::cutree(hc, k=N_clusters))
			myclusters_num <- levels(myclusters)
			if(N_clusters<9){
				levels(myclusters) <- (c(grDevices::hsv(0.56, 0.9, 1), 
										 grDevices::hsv(0, 0.27, 1), 
										 grDevices::hsv(0.52, 1, 0.5), 
										 grDevices::hsv(0.12, 0.55, 0.97), 
										 grDevices::hsv(0.83, 0.81, 0.55), 
										 grDevices::hsv(0.66, 0.15, 1),
										 grDevices::hsv(0.7, 1, 0.7), 
										 grDevices::hsv(0.42, 0.33, 1)
				)[1:N_clusters]
				)
			}else{
				levels(myclusters) <- grDevices::rainbow(N_clusters, start=0.1, end=0.9)
			}
			myclusters <- as.character(myclusters)
		}else if(!is.null(myclusters)){
			myclusters_num <- 1:length(levels(as.factor(myclusters)))
			N_clusters <- length(levels(as.factor(myclusters)))
		}
		
		if(plot){
			
			myhclustward<- function(d, method = "ward.D2", members=NULL){
				stats::hclust(d, method = "ward.D2", members=NULL)
			}
			
			if(is.null(N_clusters)){
				heatKey.size <- 2.6*heatKey.size
				dendrogram.size <- 2*dendrogram.size
			}else{
				heatKey.size <- 0.9*heatKey.size
				dendrogram.size <- 0.8*dendrogram.size
			}
			
			#d <- length(unique(map2heat))
			#if(floor(d/2)!=d/2){d <- d-1}
			maxAbs <- max(abs(min(map2heat)), abs(max(map2heat)))
			if(is.null(legend.breaks)){
				legend.breaks <- c(seq(from=-ceiling(maxAbs*10)/10, to=0, by=ceiling(maxAbs)/100),
								   seq(from=0, to=ceiling(maxAbs*10)/10, by=ceiling(maxAbs)/100)[-1]
				)
				#legend.breaks <- c(seq(from=-maxAbs, to=0, length.out=d/2+1),
				#                   seq(from=0, to=maxAbs, length.out=d/2+1)[-1]
				#)
			}
			
			colnames(map2heat) <- paste(time_unit, colnames(map2heat), sep="")
			
			if(ranking){
				clRows <- FALSE
				dendo <- "none"
			}
			else{
				dendo <- "row"
			}
			
			try(
				if(!is.null(myclusters)){
					MYheatmap.2(x = map2heat, 
								Rowv=clRows,
								Colv=FALSE,
								hclustfun=myhclustward,
								dendrogram=dendo,
								scale="none",
								col=grDevices::colorRampPalette(rev(color.vec))(length(legend.breaks)-1),
								breaks=legend.breaks,
								symkey=TRUE,
								colsep=NULL,
								rowsep=NULL, #c(1:dim(map2heat)[1]),
								sepwidth=c(0.01,0.01),
								trace="none",
								RowSideColors = myclusters,
								density.info="none",
								lmat=matrix(c(5,6,4,3,1,2), nrow=2,ncol=3,byrow=TRUE),
								lhei=c(0.115*heatKey.size, 0.3*heatmap.height),
								lwid=c(0.1*dendrogram.size,0.01,0.4*heatmap.width),
								cexRow = 0.1*cex.label.row + 0.5*cex.label.row/log10(dim(map2heat)[1]),
								cexCol = 0.1*cex.label.column + 0.5*cex.label.column/log10(dim(map2heat)[2]),
								margins=margins,
								main=list(mymain, cex=1.3*cex.main),
								...
					)
					if(N_clusters<9){
						legFill<- (c(grDevices::hsv(0.56, 0.9, 1), 
									 grDevices::hsv(0, 0.27, 1), 
									 grDevices::hsv(0.52, 1, 0.5), 
									 grDevices::hsv(0.12, 0.55, 0.97), 
									 grDevices::hsv(0.83, 0.81, 0.55), 
									 grDevices::hsv(0.66, 0.15, 1),
									 grDevices::hsv(0.7, 1, 0.7), 
									 grDevices::hsv(0.42, 0.33, 1)
						)[1:N_clusters]
						)
					}else{
						legFill <- grDevices::rainbow(N_clusters, start=0.1, end=0.9)
						legFill <- legFill[order(legFill)]
					}
					if(is.null(label.clusters) | length(label.clusters)!=N_clusters){
						graphics::legend("topright", legend=myclusters_num, fill=legFill, cex=0.7*cex.clusterKey, horiz=horiz.clusterKey)
					}else{
						graphics::legend("topright", legend=label.clusters, fill=legFill, cex=0.7*cex.clusterKey, horiz=horiz.clusterKey)
					}        
				}else{
					MYheatmap.2(x = map2heat,
								Rowv=clRows,
								Colv=FALSE,
								hclustfun=myhclustward,
								dendrogram=dendo,
								scale="none",
								col=grDevices::colorRampPalette(rev(color.vec))(length(legend.breaks)-1),
								breaks=legend.breaks,
								symkey=TRUE,
								colsep=NULL,
								rowsep=NULL, #c(1:dim(map2heat)[1]),
								sepwidth=c(0.01,0.01),
								trace="none",
								density.info="none",
								lmat=matrix(c(4,3,2,1), nrow=2,ncol=2,byrow=TRUE),
								lhei=c(0.115*heatKey.size, 0.3*heatmap.height),
								lwid=c(0.1*dendrogram.size,0.4*heatmap.width),
								cexRow = 0.1*cex.label.row + 0.5*cex.label.row/log10(dim(map2heat)[1]),
								cexCol = 0.1*cex.label.column + 0.5*cex.label.column/log10(dim(map2heat)[2]),
								margins=margins,
								main=list(mymain, cex=1.3*cex.main),
								...
					)
				}
			)
		}
		hc$legend.breaks <- legend.breaks
		hc$myclusters <- myclusters
		hc$ddr <- ddr
		hc$order <- rev(stats::order.dendrogram(ddr))
		hc$geneset.names <- gmt$geneset.names[select]
		hc$clust.trends <- clust_trends
		hc$data_plotted <- as.matrix(map2heat[hc$order,])
		
		if(!is.null(myclusters)){
			clusters <- as.factor(myclusters)
			if(!is.null(label.clusters)){
				levels(clusters) <- label.clusters
			}else{
				levels(clusters) <- myclusters_num
			}
			GSs <- rownames(map2heat)
			clustersExport <- cbind.data.frame("GeneSetTrend"=GSs, "Cluster"=as.character(clusters))
			rownames(clustersExport) <- GSs
			clustersExport  <-  clustersExport[order(clusters),]
		}else{
			clustersExport<- NULL
		}
		
		hc$clusterExport <- clustersExport
		
		invisible(hc)
	}
