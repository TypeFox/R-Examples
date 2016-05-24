#'Plotting a Specific Gene Set Stratifying on Patients
#'
#'This function can plot different representations of the gene expression in a
#'specific gene set, stratified on all subjects.
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
#'the same order as \code{TimePoint} and the columns of \code{expr} (when it is
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
#'that there is only one treatment group.  See Details.
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
#@param indiv a character string indicating by which unit observations are
#aggregated (through \code{aggreg.fun}) before the clustering.  Possible
#values are \code{"genes"} or \code{"patients"}.  Default is \code{"genes"}.
#See Details.
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
#@param showTrend 
#logical flag.  If \code{TRUE}, a black line is added for
#each cluster, representing the corresponding \code{trend.fun}.  Default is
#\code{TRUE}.
#'
#@param smooth 
#logical flag.  If \code{TRUE} and \code{showTrend} is also
#\code{TRUE}, the representation of each cluster \code{trend.fun} is smoothed
#using cubic polynoms (see \code{\link[ggplot2:stat_smooth]{stat_smooth}}.
#Default is \code{TRUE}.
#'
#'@param time_unit 
#'the time unit to be displayed (such as \code{"Y"},
#'\code{"M"}, \code{"W"}, \code{"D"}, \code{"H"}, etc) next to the values of
#'\code{TimePoint} on the x-axis.  Default is \code{""}.
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
#'@importFrom cluster agnes
#'
#'@importFrom grDevices rainbow
#'
#'@importFrom stats cutree
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
#'plotPat.1GS(expr=expr_1grp, TimePoint=design$TimePoint, 
#'        Subject_ID=design$Patient_ID, gmt=gmt_sim,
#'        geneset.name="Gene set 4",
#'        clustering=FALSE,
#'        time_unit="H",
#'        lab.cex=0.7)
#'
#'\dontrun{
#'plotPat.1GS(expr=expr_1grp, TimePoint=design$TimePoint, 
#'        Subject_ID=design$Patient_ID, gmt=gmt_sim,
#'        geneset.name="Gene set 4",
#'        clustering=FALSE, baseline=1,
#'        time_unit="H",
#'        lab.cex=0.7)
#'}        
#'
#'\dontrun{
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
#'plotPat.1GS(expr=expr_1grp, TimePoint=design$TimePoint, 
#'        Subject_ID=design$Patient_ID, gmt=gmt_sim,
#'        geneset.name="Gene set 5",
#'        time_unit="H",
#'        title="",
#'        gg.add=list(scale_color_manual(values=colval)),
#'        lab.cex=0.7
#')
#'}
#'
#'\dontrun{
#'plotPat.1GS(expr=tcgsa_sim_1grp$Estimations, TimePoint=design$TimePoint, 
#'        Subject_ID=design$Patient_ID, gmt=gmt_sim,
#'        geneset.name="Gene set 3",
#'        time_unit="H",
#'        lab.cex=0.7
#')
#'}
#'
plotPat.1GS <- 
	function(expr, gmt, Subject_ID,TimePoint, geneset.name, 
			 baseline=NULL,
			 group.var=NULL, Group_ID_paired=NULL, ref=NULL, group_of_interest=NULL,
			 FUNcluster=NULL, clustering_metric="euclidian", clustering_method="ward", B=500,
			 max_trends=4, aggreg.fun="median", trend.fun="median",
			 methodOptiClust = "firstSEmax",
			 #indiv="genes",
			 verbose=TRUE,
			 clustering=TRUE, 
			 #showTrend=TRUE, smooth=TRUE,
			 time_unit="", title=NULL, y.lab=NULL, desc=TRUE,
			 lab.cex=1, axis.cex=1, main.cex=1, y.lab.angle=90, x.axis.angle=45,
			 y.lim=NULL, x.lim=NULL, 
			 gg.add=list(theme())
	){
		
		capwords <- function(s, strict = FALSE){
			cap <- function(s){
				paste(toupper(substring(s,1,1)),{s <- substring(s,2); if(strict) tolower(s) else s},
					  sep = "", collapse = " ")
			}
			sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
		}
		
		Fun_byIndex<-function(X, index, fun){
			tapply(X, INDEX=index, FUN = fun)
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
				mytitle <- paste(geneset.name, "\n", mydesc, "\n genes \n", sep="")
				main.cex <- main.cex*0.15
				lab.cex <- lab.cex*0.5
				axis.cex <- axis.cex*0.5
			}else{
				mytitle <- paste(geneset.name, "\n genes \n", sep="")
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
		
		all_clust <- plot1GS(expr, gmt, Subject_ID, TimePoint, geneset.name, 
							 baseline, 
							 group.var, Group_ID_paired, ref, group_of_interest,
							 FUNcluster, clustering_metric, clustering_method, B,
							 max_trends, aggreg.fun, trend.fun,
							 methodOptiClust,
							 indiv="genes",
							 verbose,
							 clustering, showTrend=FALSE, smooth=FALSE, precluster=NULL,
							 time_unit, title, y.lab, desc, margins=1, line.size=1,
							 lab.cex, axis.cex, main.cex, y.lab.angle, x.axis.angle,
							 y.lim, x.lim, 
							 gg.add, 
							 plot=FALSE)
		data_stand <- t(apply(X=data_sel, MARGIN=1, FUN=scale))
		
		if(!is.null(baseline)){
			for(p in unique(Subject_ID)){
				colbaseline <- which(sort(unique(TimePoint))==baseline & Subject_ID==p)
				if(length(colbaseline)==0){
					stop("the 'baseline' value used is not one of the time points in 'TimePoint'...\n\n")
				}
				data_stand[, which(Subject_ID==p)] <- data_stand[, which(Subject_ID==p)]-data_stand[,colbaseline]
			}
		}
		
		rownames(all_clust) <- all_clust$ProbeID
		clust <- all_clust[rownames(data_stand), "Cluster"]
		
		
		meltedData <- NULL
		subj_temp <- NULL
		for (p in unique(Subject_ID)){
			sample_sel <- which(Subject_ID==p)
			melted_temp <- melt(cbind.data.frame("Probe_ID"=rownames(data_stand[, sample_sel]), 
												 "Cluster"=clust, data_stand[, sample_sel]), 
								id.vars=c("Probe_ID", "Cluster"), variable.name="TimePoint")
			meltedData <- rbind(meltedData, melted_temp)
			subj_temp <- c(subj_temp, rep(p, dim(melted_temp)[1]))
		}
		meltedData$Subject_ID <- as.factor(subj_temp)
		rm(list=c("melted_temp", "subj_temp"))
		meltedData$Cluster <- as.factor(meltedData$Cluster)
		
		
		
		
		meltedData$TimePoint <- paste(time_unit, meltedData$TimePoint, sep="")
		
		# Next TODO: medoids by patient
		# browser()
		# meltedStats$TimePoint <- paste(time_unit, meltedStats$TimePoint, sep="")
		
		if(is.null(y.lim)){
			y.max <- max(abs(meltedData$value))
			y.min <- -y.max
		}else{
			y.max <- y.lim[2]
			y.min <- y.lim[1]
		}
		if(is.null(x.lim)){
			x.lim <- unique(meltedData$TimePoint)
		}
		p <- (ggplot(meltedData, aes_string(x="TimePoint", y="value")) 
			  + geom_hline(aes(yintercept = 0), linetype=1, colour='grey50', size=0.4)
			  + facet_wrap( ~Subject_ID, ncol=floor(sqrt(length(unique(meltedData$Subject_ID)))))
		)
		
		if(!clustering){
			p <- (p
				  + geom_line(aes_string(group="Probe_ID", colour="Probe_ID"), size=0.7)
				  + scale_colour_manual(guide='none', name='probe ID', 
				  					  values=grDevices::rainbow(length(select_probe)))
			)
		}else{
			p <- (p
				  + geom_line(aes_string(group="Probe_ID", colour="Cluster"), size=0.7)
				  + guides(colour = guide_legend(override.aes=list(size=1, fill="white"), keywidth=2*lab.cex, 
				  							   title.theme=element_text(size = 15*lab.cex, angle=0),
				  							   label.theme=element_text(size = 9*lab.cex, angle=0)
				  )
				  )
			)
		}
		
		p <- (p
			  + ylab(y.lab)
			  + xlab('Time')
			  + ggtitle(mytitle)
			  + theme(plot.title=element_text(size = 35*main.cex))
			  + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_rect(colour='grey40', fill = 'white'))
			  + ylim(y.min, y.max)
			  + xlim(x.lim)
			  + theme(axis.title.y = element_text(size = 25*lab.cex, angle = y.lab.angle, vjust=0.3), axis.text.y = element_text(size=18*axis.cex, colour = 'grey40')) 
			  + theme(axis.title.x = element_text(size = 25*lab.cex, angle = 0, vjust=-0.9), axis.text.x = element_text(size=18*axis.cex, colour = 'grey40', angle=x.axis.angle, vjust=0.5, hjust=0.5))
			  + theme(plot.margin=unit(c(0.5, 0.5, 0.7, 1), 'lines'))
			  + theme(legend.key=element_rect(fill="white"))
		)
		
		for(a in gg.add){
			p <- p + a
		}
		
		#   if(showTrend){
		#     if(!smooth){
		#       p <- (p + geom_line(data=meltedStats, aes(x=TimePoint, y=value, group=Cluster, linetype=Cluster), size=4))
		#     }else{
		#       p <- (p + stat_smooth(formula=y~poly(x,3), data=meltedStats, aes(x=TimePoint, y=value, group=Cluster, linetype=Cluster), size=4, se=FALSE, method="lm", color="black"))
		#     }
		#     p <- p + scale_linetype_manual(name=paste("Cluster", capwords(trend.fun)), values=as.numeric(levels(meltedStats$Cluster))+1, 
		#                                      guide=guide_legend(override.aes=list(size=1), keywidth=2*lab.cex, 
		#                                                         title.theme=element_text(size = 17*lab.cex, angle=0),
		#                                                         label.theme=element_text(size = 12*lab.cex, angle=0)
		#                                      )
		#     )
		#   }
		print(p)
		invisible(clust)
	}



