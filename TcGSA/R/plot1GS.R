#'Plotting a Specific Gene Set
#'
#'This function can plot different representations of the gene expression in a
#'specific gene set.
#'
#'If \code{expr} is a matrix or a dataframe, then the "original" data are
#'plotted.  On the other hand, if \code{expr} is a list returned in the
#'\code{'Estimations'} element of \code{\link{TcGSA.LR}}, then it is those
#'"estimations" made by the \code{\link{TcGSA.LR}} function that are plotted.
#'
#'If \code{indiv} is 'genes', then each line of the plot is the median of a
#'gene expression over the patients. On the other hand, if \code{indiv} is
#''patients', then each line of the plot is the median of a patient genes
#'expression in this gene set.
#'
#'This function uses the Gap statistics to determine the optimal number of
#'clusters in the plotted gene set.  See
#'\code{\link[cluster:clusGap]{clusGap}}.
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
#'@param gmt 
#'a \bold{gmt} object containing the gene sets definition.  See
#'\code{\link[GSA:GSA.read.gmt]{GSA.read.gmt}} and
#'definition on \href{http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats}{www.broadinstitute.org}.
#'
#'@param Subject_ID 
#'a factor of length \eqn{p} that is in the same order as the
#'columns of \code{expr} (when it is a dataframe) and that contains the patient
#'identifier of each sample.
#TODO See Details.
#'
#'@param TimePoint 
#'a numeric vector or a factor of length \eqn{p} that is in
#'the same order as \code{Subject_ID} and the columns of \code{expr} (when it is
#'a dataframe), and that contains the time points at which gene expression was
#'measured.
#TODO See Details.
#'
#'@param geneset.name 
#'a character string containing the name of the gene set to
#'be plotted, that must appear in the \code{"geneset.names"} element of
#'\code{gmt}.
#'
#'@param baseline 
#'a character string which is the value of \code{TimePoint}
#'that can be used as a baseline.  Default is \code{NULL}, in which case no
#'timepoint is used as a baseline value for gene expression.  Has to be
#'\code{NULL} when comparing two treatment groups.
#TODO See Details.
#'
#'@param group.var 
#'in the case of several treatment groups, this is a factor of
#'length \eqn{p} that is in the same order as \code{Timepoint},
#'\code{Subject_ID} and the columns of \code{expr}.  It indicates to which
#'treatment group each sample belongs to.  Default is \code{NULL}, which means
#'that there is only one treatment group.  
#TODO See Details.
#'
#'@param Group_ID_paired 
#'a character vector of length \eqn{p} that is in the
#'same order as \code{Timepoint}, \code{Subject_ID}, \code{group.var} and the
#'columns of \code{expr}.  This argument must not be \code{NULL} in the case of
#'a paired analysis, and must be \code{NULL} otherwise.  Default is
#'\code{NULL}.  
#TODO See Details.
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
#'of the labels of \code{group.var}.  
#TODO See Details.
#'
#'@param FUNcluster 
#'a function which accepts as first argument a matrix
#'\code{x} and as second argument the number of clusters desired \code{k}, and
#'which returns a list with a component named \code{'cluster'} which is a
#'vector of length n = nrow(x) of integers in 1:k, determining the clustering
#'or grouping of the n observations.  Default is \code{NULL}, in which case a
#'hierachical clustering is performed via the function
#'\code{\link[cluster:agnes]{agnes}}, using the metric \code{clustering_metric}
#'and the method \code{clustering_method}.  See \code{'FUNcluster'} in
#'\code{\link[cluster:clusGap]{clusGap}} and Details.
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
#'or the name of any other defined statistics function that returns a single
#'numeric value.  It specifies the function used to aggregate the observations
#'before the clustering.  Default is to \code{median}.
#'
#'@param trend.fun 
#'a character string such as \code{"mean"}, \code{"median"} or
#'the name of any other function that returns a single numeric value.  It
#'specifies the function used to calculate the trends of the identified
#'clustered.  Default is to \code{median}.
#'
#'@param methodOptiClust 
#'character string indicating how the "optimal" number
#'of clusters is computed from the gap statistics and their standard
#'deviations. Possible values are \code{"globalmax"}, \code{"firstmax"},
#'\code{"Tibs2001SEmax"}, \code{"firstSEmax"} and \code{"globalSEmax"}.
#'Default is \code{"firstSEmax"}.  See \code{'method'} in
#'\code{\link[cluster:clusGap]{clusGap}}, Details and \emph{Tibshirani et al.,
#'2001} in References.
#'
#'@param indiv a character string indicating by which unit observations are
#'aggregated (through \code{aggreg.fun}) before the clustering.  Possible
#'values are \code{"genes"} or \code{"patients"}.  Default is \code{"genes"}.
#'See Details.
#'
#'@param verbose 
#'logical flag enabling verbose messages to track the computing
#'status of the function.  Default is \code{TRUE}.
#'
#'@param clustering 
#'logical flag.  If \code{FALSE}, there is no clustering
#'representation; if \code{TRUE}, the lines are colored according to which
#'cluster they belong to.  Default is \code{TRUE}.  See Details.
#'
#'@param showTrend 
#'logical flag.  If \code{TRUE}, a black line is added for
#'each cluster, representing the corresponding \code{trend.fun}.  Default is
#'\code{TRUE}.
#'
#'@param smooth 
#'logical flag.  If \code{TRUE} and \code{showTrend} is also
#'\code{TRUE}, the representation of each cluster \code{trend.fun} is smoothed
#'using cubic polynoms (see \code{\link[ggplot2:geom_smooth]{geom_smooth}}.
#'Default is \code{TRUE}. 
#'At the moment, must accept parameter \code{"na.rm"} (which is automatically set to \code{TRUE}). 
#'This might change in future versions
#'
#'@param precluster 
#'a vector of length \eqn{p} that is in
#'the same order as \code{Subject_ID}, \code{TimePoint} and the columns of \code{expr} (when it is
#'a dataframe), and that contains a prior clustering of the subjects. Default is \code{NULL}.
#'
#'@param time_unit 
#'the time unit to be displayed (such as \code{"Y"},
#'\code{"M"}, \code{"W"}, \code{"D"}, \code{"H"}, etc) next to the values of
#'\code{TimePoint} on the x-axis.  Default is \code{""}, in which case the time 
#'scale on the x-axis is proportionnal to the time values.
#'
#'@param title 
#'character specifying the title of the plot.  If \code{NULL}, a
#'title is automatically generated, if \code{""}, no title appears.  Default is
#'\code{NULL}.
#'
#'@param y.lab 
#'character specifying the annotation of the y axis.  If \code{NULL}, an
#'annotation is automatically generated, if \code{""}, no annotation appears.  Default is
#'\code{NULL}.
#'
#'@param desc 
#'a logical flag. If \code{TRUE}, a line is added to the title of
#'the plot with the description of the gene set plotted (from the gmt file).
#'Default is \code{TRUE}.
#'
#'@param lab.cex 
#'a numerical value giving the amount by which lab labels text
#'should be magnified relative to the default \code{1}.
#'
#'@param axis.cex 
#'a numerical value giving the amount by which axis annotation
#'text should be magnified relative to the default \code{1}.
#'
#'@param main.cex 
#'a numerical value giving the amount by which title text
#'should be magnified relative to the default \code{1}.
#'
#'@param margins
#'a numerical value giving the amount by which the margins
#'should be reduced or increased relative to the default \code{1}.
#'
#'@param line.size
#'a numerical value giving the amount by which the line sizes
#'should be reduced or increased relative to the default \code{1}.
#'
#'@param y.lab.angle 
#'a numerical value (in [0, 360]) giving the orientation by
#'which y-label text should be turned (anti-clockwise).  Default is \code{90}.
#'See \code{\link[ggplot2:element_text]{element_text}}.
#'
#'@param x.axis.angle 
#'a numerical value (in [0, 360]) giving the orientation by
#'which x-axis annotation text should be turned (anti-clockwise).  Default is
#'\code{45}.
#'
#'@param y.lim 
#'a numeric vector of length 2 giving the range of the y-axis.
#'See \code{\link{plot.default}}.
#'
#'@param x.lim 
#'if numeric, will create a continuous scale, if factor or
#'character, will create a discrete scale.  Observations not in this range will
#'be dropped.  See \code{\link{xlim}}.
#'
#'@param gg.add 
#'A list of instructions to add to the ggplot2 instruction.  See \link{+.gg}.  Default is \code{list(theme())}, which adds nothing
#'to the plot.
#'
#'@param plot 
#'logical flag.  If \code{FALSE}, no plot is drawn.  Default is \code{TRUE}.
#'
#'@return A dataframe the 2 following variables: \itemize{
#'\item \code{ProbeID} which contains the IDs of the probes of the plotted gene set.
#'\item \code{Cluster} which to which cluster the probe belongs to.
#'} If \code{clustering} is \code{FALSE}, then \code{Cluster} is \code{NA} for all the probes.
#'
#'@author Boris P. Hejblum
#'
#'@seealso \code{\link[ggplot2:ggplot]{ggplot}}, \code{\link[cluster:clusGap]{clusGap}}
#'
#'@references Tibshirani, R., Walther, G. and Hastie, T., 2001, Estimating the
#'number of data clusters via the Gap statistic, \emph{Journal of the Royal
#'Statistical Society, Series B (Statistical Methodology)}, \bold{63}, 2:
#'411--423.
#'
#'@import ggplot2
#'
#'@importFrom cluster agnes clusGap maxSE
#'
#'@importFrom stats cutree lm
#'
#'@export
#'
#'@examples
#'
#'data(data_simu_TcGSA)
#'tcgsa_sim_1grp <- TcGSA.LR(expr=expr_1grp, gmt=gmt_sim, design=design, 
#'                           subject_name="Patient_ID", time_name="TimePoint",
#'                           time_func="linear", crossedRandom=FALSE)
#'
#'plot1GS(expr=expr_1grp, TimePoint=design$TimePoint, 
#'        Subject_ID=design$Patient_ID, gmt=gmt_sim,
#'        geneset.name="Gene set 4",
#'        indiv="genes", clustering=FALSE,
#'        time_unit="H",
#'        lab.cex=0.7)
#'
#'\dontrun{ 
#'plot1GS(expr=expr_1grp, TimePoint=design$TimePoint, 
#'        Subject_ID=design$Patient_ID, gmt=gmt_sim,
#'        geneset.name="Gene set 5",
#'        indiv="patients", clustering=FALSE, baseline=1,
#'        time_unit="H",
#'        lab.cex=0.7)
#'}
#'\dontrun{        
#'plot1GS(expr=tcgsa_sim_1grp$Estimations, TimePoint=design$TimePoint, 
#'        Subject_ID=design$Patient_ID, gmt=gmt_sim,
#'        geneset.name="Gene set 5",
#'        indiv="genes",
#'        time_unit="H",
#'        lab.cex=0.7
#')
#'}
#'
#'\dontrun{
#'library(grDevices)
#'library(graphics)
#'colval <- c(hsv(0.56, 0.9, 1),
#'            hsv(0, 0.27, 1),
#'            hsv(0.52, 1, 0.5),
#'            hsv(0, 0.55, 0.97),
#'            hsv(0.66, 0.15, 1),
#'            hsv(0, 0.81, 0.55),
#'            hsv(0.7, 1, 0.7),
#'            hsv(0.42, 0.33, 1)
#')
#'n <- length(colval);  y <- 1:n
#'op <- par(mar=rep(1.5,4))
#'plot(y, axes = FALSE, frame.plot = TRUE,
#'		 xlab = "", ylab = "", pch = 21, cex = 8,
#'		 bg = colval, ylim=c(-1,n+1), xlim=c(-1,n+1),
#'		 main = "Color scale"
#')
#'par(op)
#'
#'require(ggplot2)
#'plot1GS(expr=expr_1grp, TimePoint=design$TimePoint, 
#'        Subject_ID=design$Patient_ID, gmt=gmt_sim,
#'        geneset.name="Gene set 5",
#'        indiv="genes",
#'        time_unit="H",
#'        title="",
#'        gg.add=list(scale_color_manual(values=colval), 
#'                    guides(colour = guide_legend(reverse=TRUE))),
#'        lab.cex=0.7
#')
#'}
#'
plot1GS <- 
	function(expr, gmt, Subject_ID, TimePoint, geneset.name, 
			 baseline=NULL,
			 group.var=NULL, Group_ID_paired=NULL, ref=NULL, group_of_interest=NULL,
			 FUNcluster=NULL, clustering_metric="euclidian", clustering_method="ward", B=500,
			 max_trends=4, aggreg.fun="median", trend.fun="median",
			 methodOptiClust = "firstSEmax",
			 indiv="genes",
			 verbose=TRUE,
			 clustering=TRUE, showTrend=TRUE, smooth=TRUE, precluster=NULL, 
			 time_unit="", title=NULL, y.lab=NULL, desc=TRUE,
			 lab.cex=1, axis.cex=1, main.cex=1, y.lab.angle=90, x.axis.angle=45, margins=1, line.size=1,
			 y.lim=NULL, x.lim=NULL, 
			 gg.add=list(theme()),
			 plot=TRUE
	){
		
		pre_clustering <- !is.null(precluster)
		
		capwords <- function(s, strict = FALSE){
			cap <- function(s){
				paste(toupper(substring(s,1,1)),{s <- substring(s,2); if(strict) tolower(s) else s},
					  sep = "", collapse = " ")
			}
			sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
		}
		
		Fun_byIndex<-function(X, index, fun, ...){
			tapply(X, INDEX=index, FUN = fun, ...)
		}
		
		if(is.null(FUNcluster)){
			FUNcluster <- switch(EXPR=clustering_metric,
								 sts= function(x, k, time, ...){
								 	d <- STSdist(m=x, time = time)
								 	clus <- stats::cutree(agnes(d, ...), k=k)
								 	return(list("cluster"=clus))
								 },
								 function(x, k, ...){
								 	clus <- stats::cutree(agnes(x, method=clustering_method, metric=clustering_metric, ...), k=k)
								 	return(list("cluster"=clus))
								 }
			)
		}
		if(!is.function(FUNcluster)){
			stop("the 'FUNcluster' supplied is not a function")
		}
		
		if(is.null(title)){
			if(desc){
				mydesc <- gmt$geneset.descriptions[which(gmt$geneset.names==geneset.name)]
				mytitle <- paste(geneset.name, "\n", mydesc, "\n", indiv, "\n", sep="")
				main.cex <- main.cex*0.15
				lab.cex <- lab.cex*0.5
				axis.cex <- axis.cex*0.5
			}else{
				mytitle <- paste(geneset.name, "\n", indiv, "\n", sep="")
			}
		}else{
			if(title==""){
				mytitle <- NULL
			}else{
				mytitle <- title
			}
		}
		
		if(is.null(y.lab)){
			if(is.data.frame(expr)){
				y.lab <- paste(capwords(aggreg.fun), 'of standardized gene expression')
			}
			else{
				y.lab <- paste(capwords(aggreg.fun), 'of standardized estimate')
			}
		}
		
		
		interest <- which(gmt$geneset.names==geneset.name)
		if(length(interest)==0){
			stop("The 'geneset.name' supplied is not in the 'gmt'")
		}
		if(is.data.frame(expr)){
			select_probe <- intersect(rownames(expr), unique(gmt$genesets[[interest]]))
			data_sel <- as.matrix(expr[select_probe,])
		}else if(is.list(expr)){
			expr_sel <- expr[[interest]]
			expr_sel <- expr_sel[, , order(as.numeric(dimnames(expr_sel)[[3]])), drop=FALSE]
			data_sel <- matrix(expr_sel, nrow=dim(expr_sel)[1], ncol=dim(expr_sel)[2]*dim(expr_sel)[3])
			rownames(data_sel) <- dimnames(expr_sel)[[1]]
			select_probe <- dimnames(expr_sel)[[1]]
			TimePoint <- sort(as.numeric(rep(dimnames(expr_sel)[[3]], dim(expr_sel)[2])))
			Subject_ID <- rep(dimnames(expr_sel)[[2]], dim(expr_sel)[3])
		}
		
		data_stand <- t(apply(X=data_sel, MARGIN=1, FUN=scale))
		if(indiv=="genes"){
			data_stand_MedianByTP <- t(apply(X=data_stand, MARGIN=1, FUN=Fun_byIndex, index=as.factor(TimePoint), fun=aggreg.fun, na.rm=T))
		}else if(indiv=="patients"){
			data_tocast<-cbind.data.frame(TimePoint, Subject_ID, "M" = apply(X=data_stand, MARGIN=2, FUN=aggreg.fun))
			data_stand_MedianByTP <- as.matrix(acast(data_tocast, formula="Subject_ID~TimePoint", value.var="M"))
		}
		
		
		if(!is.null(baseline)){
			colbaseline <- which(sort(unique(TimePoint))==baseline)
			if(length(colbaseline)==0){
				stop("the 'baseline' value used is not one of the time points in 'TimePoint'...\n\n")
			}
			data_stand_MedianByTP <- data_stand_MedianByTP-data_stand_MedianByTP[,colbaseline]
		}
		
		if(!pre_clustering && (clustering | showTrend) && length(which(is.na(data_stand_MedianByTP)))>0){
			warning("Unable to compute optimal number of clusters/trends due to missing values\n")
			clustering <- FALSE
			showTrend <- FALSE
			
		}
		
		if(clustering | showTrend){
			if(!pre_clustering){
				if(verbose){
					cat("Optimally clustering...\n")
				}
				kmax <- ifelse(dim(data_stand_MedianByTP)[1]>4, max_trends, dim(data_stand_MedianByTP)[1]-1)
				
				if(kmax>=2){
					if(clustering_metric!="sts"){
						cG <- clusGap(x=data_stand_MedianByTP, FUNcluster=FUNcluster, K.max=kmax, B=B, verbose=FALSE)
						nc <- maxSE(f = cG$Tab[, "gap"], SE.f = cG$Tab[, "SE.sim"], method = methodOptiClust)
						clust <- FUNcluster(data_stand_MedianByTP, k=nc)$cluster
					}else{
						cG <- clusGap(x=data_stand_MedianByTP, FUNcluster=FUNcluster, K.max=kmax, B=B, verbose=FALSE, time=as.numeric(colnames(data_stand_MedianByTP)))
						nc <- maxSE(f = cG$Tab[, "gap"], SE.f = cG$Tab[, "SE.sim"], method = methodOptiClust)
						clust <- FUNcluster(data_stand_MedianByTP, k=nc, time=as.numeric(colnames(data_stand_MedianByTP)))$cluster
					}    
				}else{
					nc <- 1
					clust <- rep(1, dim(data_stand_MedianByTP)[1])
				}
				
				
				
				medoids <- as.data.frame(t(apply(X=data_stand_MedianByTP, MARGIN=2, FUN=Fun_byIndex, index=clust, fun=trend.fun)))
				if(dim(medoids)[1]==1){
					medoids <- cbind.data.frame("TimePoint"= colnames(medoids), "1"=t(medoids))
				}else{
					medoids <- cbind.data.frame("TimePoint"= rownames(medoids), medoids)
				}
			}else{
				
				clust <- precluster
				medoids <- as.data.frame(t(apply(X=data_stand_MedianByTP, MARGIN=2, FUN=Fun_byIndex, 
												 index=as.factor(as.numeric(precluster)), fun=trend.fun, na.rm=TRUE)))
				if(nrow(medoids)==1){
					medoids <- cbind.data.frame("TimePoint"= colnames(medoids), "1"=t(medoids))
				}else{
					medoids <- cbind.data.frame("TimePoint"= rownames(medoids), medoids)
				}
				colnames(medoids) <- c("TimePoint", levels(precluster))
			}
			if(verbose){
				cat("DONE\n")
			}
			
			
			if(indiv=="patients"){
				row_names <- (as.numeric(sub("P","",rownames(data_stand_MedianByTP))))
				position<-NULL
				for(i in 1:length(row_names)){
					for(j in 1:length(row_names)){
						if(row_names[i]==row_names[j]){
							position<-c(position,j)
						}
					}
				}
				position<-order(as.character(position))
			}
		}else{
			medoids <- cbind.data.frame("TimePoint"=colnames(data_stand_MedianByTP), "1"='NA')
			clust <- rep(NA, dim(data_stand_MedianByTP)[1])
		}
		
		medoids$TimePoint <- as.numeric(as.character(medoids$TimePoint))
		colnames(data_stand_MedianByTP) <- as.numeric(colnames(data_stand_MedianByTP))
		
		classif <- cbind.data.frame("ProbeID"=rownames(data_stand_MedianByTP), "Cluster"=clust)
		
		meltedData <- melt(cbind.data.frame("Probe_ID"=rownames(data_stand_MedianByTP), "Cluster"=classif$Cluster, data_stand_MedianByTP), id.vars=c("Probe_ID", "Cluster"), variable.name="TimePoint")
		meltedStats <- melt(medoids, id.vars="TimePoint", variable.name="Cluster")
		meltedData$Cluster <- as.character(meltedData$Cluster)
		meltedData$TimePoint <- as.numeric(as.character(meltedData$TimePoint))
		meltedStats$TimePoint <- as.numeric(as.character(meltedStats$TimePoint))
		MeasPt <- unique(meltedData$TimePoint)
		
		
		# 		browser()
		# 		pca_data <- acast(meltedData, formula=TimePoint~Probe_ID)
		# 		pca_res <- dudi.pca(pca_data)
		# 		1
		# 		pdf(width=5, height=4, file="~/PCAc1_M6_7.pdf")
		# 		plot(y=-pca_res$l1[,1], x=rownames(pca_res$l1), type="l",
		# 			 ylab="PCA 1st comp", xlab="Time (weeks)", main=g,
		# 			 lwd=3, col="blue")
		# 		dev.off()
		
		
		
		if(is.null(y.lim)){
			y.max <- max(abs(meltedData$value), na.rm = TRUE)
			y.min <- -y.max
		}else{
			y.max <- y.lim[2]
			y.min <- y.lim[1]
		}
		if(is.null(x.lim)){
			x.lim <- c(min(MeasPt), max(MeasPt))
		}
		
		p <- (ggplot(meltedData, aes_string(x="TimePoint", y="value")) 
			  + geom_hline(aes(yintercept = 0), linetype=1, colour='grey50', size=0.4*line.size)
			  + theme(panel.border=element_rect(fill=NA, size=0.1*line.size, colour='grey50'),
			  		axis.ticks=element_line(size=0.4*line.size, colour='grey50'))
		)
		
		myalpha <- 1 
		if(showTrend){
			myalpha <- 0.5
		}
		
		if(clustering | pre_clustering){
			p <- (p
				  + geom_line(aes_string(group="Probe_ID", colour="Cluster"), size=0.5*line.size, alpha=myalpha)
				  + guides(colour = guide_legend(override.aes=list(size=1, fill="white"), keywidth=2*lab.cex, 
				  							   title.theme=element_text(size = 15*lab.cex, angle=0),
				  							   label.theme=element_text(size = 9*lab.cex, angle=0)
				  )
				  )
			)
		}else{
			p <- (p
				  + geom_line(aes_string(group="Probe_ID", colour="Probe_ID"), size=0.5*line.size, alpha=myalpha)
			)
			if(indiv=="patients"){
				p <- (p + guides(colour=guide_legend(title='Subject')))
				#+ scale_colour_manual(guide='none', name='Subject', values=rainbow(length(select_probe))))
			}
		}
		
		p <- (p
			  + ylab(y.lab)
			  + xlab('Time')
			  + ggtitle(mytitle)
			  + theme(plot.title=element_text(size = 35*main.cex))
			  + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_rect(colour='grey40', fill = 'white'))
			  + ylim(y.min, y.max)
			  + scale_x_continuous(breaks=MeasPt, 
			  					 labels=paste(time_unit, MeasPt, sep="")
			  )
			  + theme(axis.title.y = element_text(size = 25*lab.cex, angle = y.lab.angle, vjust=0.8), axis.text.y = element_text(size=18*axis.cex, colour = 'grey40')) 
			  + theme(axis.title.x = element_text(size = 25*lab.cex, angle = 0, vjust=0.5), axis.text.x = element_text(size=18*axis.cex, colour = 'grey40', angle=x.axis.angle, vjust=0.5, hjust=0.5))
			  + theme(plot.margin=unit(margins*c(0.5, 0.7, 0.1, 0.5), 'lines'))
			  + theme(legend.key=element_rect(fill="white"))
		)
		
		
		if(showTrend){
			
			if(!smooth){
				if(pre_clustering){
					p <- (p + geom_line(data=meltedStats, aes_string(x="TimePoint", y="value", group="Cluster", colour="Cluster"), size=3))
				}else{
					p <- (p + geom_line(data=meltedStats, aes_string(x="TimePoint", y="value", group="Cluster", linetype="Cluster"), colour="black", size=3))
				}
			}else{
				if(length(MeasPt)<4){
					stop("Not enough time points to estimate a smoothed trend! 
    				 Set 'smooth' argument to 'FALSE'.\n")
				}
				# spline_knots <- ifelse(length(MeasPt)>6, 
				# 						floor((length(MeasPt)-2)/2), 
				# 						1
				# )
				if(pre_clustering){
					p <- (p 
						  + geom_smooth(formula=y~poly(x, 3), data=meltedStats, aes_string(x="TimePoint", y="value", group="Cluster", colour="Cluster", size="1.5"), se=FALSE, method="lm")
						  # + geom_smooth(data=meltedStats, aes_string(x="TimePoint", y="value", group="Cluster", colour="Cluster"), size=1.7*line.size, se=FALSE, method="loess")
						  + guides(size="none")
					)
				}else{
					p <- p + geom_smooth(formula=y~poly(x, 3), data=meltedStats, aes_string(x="TimePoint", y="value", group="Cluster", linetype="Cluster"), size=1.7*line.size, colour="black", se=FALSE, method="lm")
						  # + geom_smooth(data=meltedStats, aes_string(x="TimePoint", y="value", group="Cluster", linetype="Cluster"), size=1.7*line.size, colour="black", se=FALSE, method="loess")
					if(!clustering){
						p <- p + scale_linetype_manual(name=paste("Cluster", capwords(trend.fun)), values=as.numeric(levels(meltedStats$Cluster))+1)
						
					}
					#y~ns(x, knots=spline_knots)
				}
				if(length(unique(meltedStats$Cluster))==1 | length(unique(meltedStats$Group))==1){
					p <- (p + guides(linetype="none"))
				}
				
				if(clustering){
					p <- (p + guides(size="none")
						  + scale_linetype_manual(name=paste("Cluster", capwords(trend.fun)), values=as.numeric(levels(meltedStats$Cluster))+1,#as.numeric(levels(meltedStats$Cluster))+1, 
						  						guide=guide_legend(override.aes=list(size=1), keywidth=2*lab.cex, 
						  										   title.theme=element_text(size = 17*lab.cex, angle=0),
						  										   label.theme=element_text(size = 12*lab.cex, angle=0)
						  						)
						  )
					)
				}else{
					p <- p + scale_size_continuous(name="", labels=capwords(trend.fun))
				}	
			}
		}
		
		for(a in gg.add){
			p <- p + a
		}
		if(plot){
			print(p)
		}
		
		classif <- classif[order(classif$Cluster), ]
		invisible(classif)
		
	}



