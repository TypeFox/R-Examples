#'Cluster the genes dynamics into different dominant trends.
#'
#'This function clusters the genes dynamics of one gene sets into different
#'dominant trends.  The optimal number of custers is computed thanks to the gap
#'statistics.  See \code{\link[cluster:clusGap]{clusGap}}.
#'
#'If \code{expr} is a matrix or a dataframe, then the genes dynamics are
#'clustered on the "original" data.  On the other hand, if \code{expr} is a
#'list returned in the \code{'Estimations'} element of \code{\link{TcGSA.LR}},
#'then the dynamics are computed on the estimations made by the
#'\code{\link{TcGSA.LR}} function.
#'
#'This function uses the Gap statistics to determine the optimal number of
#'clusters in the plotted gene set.  See
#'\code{\link[cluster:clusGap]{clusGap}}.
#'
#'@aliases clustTrend ClusteredTrends print.ClusteredTrends
#'plot.ClusteredTrends
#'
#'@param tcgs a \bold{tcgsa} object for \code{clustTrend}, or a
#'\bold{ClusteredTrends} object for \code{print.ClusteredTrends} and
#'\code{plot.ClusteredTrends}.
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
#'identifier of each sample. 
#TODO See Details.
#'
#'@param TimePoint 
#'a numeric vector or a factor of length \eqn{p} that is in
#'the same order as \code{Subject_ID} and the columns of \code{expr} (when it
#'is a dataframe), and that contains the time points at which gene expression
#'was measured. 
#TODO See Details.
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
#'"\code{BH}", "\code{BY}", "\code{ABH}", "\code{TSBH}". See
#'\code{\link[multtest:mt.rawp2adjp]{mt.rawp2adjp}} for details.  Default is
#'"\code{BY}", the Benjamini & Yekutieli (2001) step-up FDR-controlling
#'procedure (general dependency structures).  In order to control the FWER(in
#'case of an analysis that is more a hypothesis confirmation than an
#'exploration of the expression data), we recommand to use "\code{Holm}", the
#'Holm (1979) step-down adjusted p-values for strong control of the FWER.
#'
#'@param nbsimu_pval 
#'the number of observations under the null distribution to
#'be generated in order to compute the p-values.  Default is \code{1e+06}.
#'
#'@param baseline 
#'a character string which is the value of \code{TimePoint}
#'that can be used as a baseline.  Default is \code{NULL}, in which case no
#'timepoint is used as a baseline value for gene expression.  Has to be
#'\code{NULL} when comparing two treatment groups.
#TODO See Details.
#'
#'@param only.signif 
#'logical flag for analysing the trends in only the
#'significant gene sets.  If \code{FALSE}, all the gene sets from the
#'\bold{gmt} object contained in \code{x} are clustered.  Default is
#'\code{TRUE}.
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
#'same order as \code{Timepoint},  \code{Subject_ID}, \code{group.var} and the
#'columns of \code{expr}.  This argument must not be \code{NULL} in the case of
#'a paired analysis, and must be \code{NULL} otherwise.  Default is
#'\code{NULL}.
#TODO See Details.
#'
#'@param ref 
#'the group which is used as reference in the case of several
#'treatment groups.  Default is \code{NULL}, which means that reference is the
#'first group in alphabetical order of the labels of \code{group.var}.
#TODO See Details.
#'
#'@param group_of_interest 
#'the group of interest, for which dynamics are to be
#'computed in the case of several treatment groups.  Default is \code{NULL},
#'which means that group of interest is the second group in alphabetical order
#'of the labels of \code{group.var}.
#TODO See Details.
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
#'or the name of any other defined statistics function that returns a single
#'numeric value.  It specifies the function used to aggregate the observations
#'before the clustering.  Default is to \code{median}.  Default is
#'\code{"median"}.
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
#'@param indiv 
#'a character string indicating by which unit observations are
#'aggregated (through \code{aggreg.fun}) before the clustering.  Possible
#'values are \code{"genes"} or \code{"patients"}.  Default is \code{"genes"}.
#'
#'@param verbose 
#'logical flag enabling verbose messages to track the computing
#'status of the function.  Default is \code{TRUE}.
#'
#'@return An object of class \bold{\link{ClusteredTrends}} which is a list with
#'the 4 following components: \itemize{
#'\item NbClust a vector that contains the optimal number of clusters for
#'each analysed gene sets.
#'\item ClustsMeds a list of the same length as \code{NsClust} (the
#'number of analysed gene sets). Each element of the list is a data frame, in
#'which there is as many column as the optimal number of clusters for the
#'corresponding gene setsfor each cluster.  Each column of the data frame
#'contains the median trend values for the corresponding cluster.
#'\item GenesPartition a list of the same length as \code{NsClust} (the
#'number of analysed gene sets).  Each element of the list is a vector which
#'gives the partition of the genes inside the corresponding gene set.
#'\item MaxNbClust an integer storing the maximum number of different
#'clusters tested, as given by the argument \code{'max_trends'}.
#'}
#'
#'@author Boris P. Hejblum
#'
#'@seealso \code{\link{plot1GS}}, \code{\link{TcGSA.LR}},
#'\code{\link[cluster:clusGap]{clusGap}}
#'
#'@references Tibshirani, R., Walther, G. and Hastie, T., 2001, Estimating the
#'number of data clusters via the Gap statistic, \emph{Journal of the Royal
#'Statistical Society, Series B (Statistical Methodology)}, \bold{63}, 2:
#'41--423.
#'
#'@importFrom cluster agnes clusGap maxSE
#'
#'@importFrom graphics barplot
#'
#'@importFrom grDevices rainbow
#'
#'@importFrom stats cutree var
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
#'  
#'CT <- clustTrend(tcgsa_sim_1grp,
#'     expr=expr_1grp, Subject_ID=design$Subject_ID, TimePoint=design$TimePoint)
#'CT
#'plot(CT)
#'
#'CT$NbClust
#'CT$NbClust["Gene set 5"]
#'CT$ClustMeds[["Gene set 4"]]
#'CT$ClustMeds[["Gene set 5"]]
#'




clustTrend <- 
	function(tcgs,
			 expr, Subject_ID, TimePoint, threshold = 0.05, 
			 myproc = "BY", nbsimu_pval = 1e+06, baseline=NULL, 
			 only.signif=TRUE, group.var=NULL, Group_ID_paired=NULL, 
			 ref=NULL, group_of_interest=NULL, FUNcluster=NULL, 
			 clustering_metric="euclidian", clustering_method="ward", 
			 B=100,	 max_trends=4, aggreg.fun="median", 
			 trend.fun="median", methodOptiClust = "firstSEmax",
			 indiv="genes", verbose=TRUE
	){
		
		
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
			# Kmeans: (be careful about about missing data)
			#     FUNcluster <- function(x, k, ...){
			#     	clus <- kmeans(x, centers=k, nstart=1)$cluster
			#     	return(list("cluster"=clus))
			#     }
		}
		if(!is.function(FUNcluster)){
			stop("the 'FUNcluster' supplied is not a function")
		}
		
		gmt <- tcgs[["GeneSets_gmt"]]
		separateSubjects <- tcgs[["separateSubjects"]]
		if(only.signif){
			GSsig <- signifLRT.TcGSA(tcgsa=tcgs, threshold=threshold, myproc = "BY", nbsimu_pval = 1e+06)$mixedLRTadjRes
			GeneSetsList <- GSsig$GeneSet
			if(length(GeneSetsList)<1){
				stop("NO SIGNIFICANT GENE SETS\n No gene sets to be plotted: set 'only.signif' argument to 'FALSE' in order to plot all the investigated gene sets")
			}
		}
		else{
			GeneSetsList <- gmt$geneset.names
		}
		
		if(!is.null(group.var)){
			if(is.null(ref)){
				ref <- levels(group.var)[1]
			}
			if(is.null(group_of_interest)){
				group_of_interest <- levels(group.var)[2]
			}
		}
		
		NbClust <- numeric(length(GeneSetsList))
		names(NbClust) <- GeneSetsList
		ClustsMeds <- list()
		GenesPartition <- list()
		for(gs in GeneSetsList){
			interest <- which(gmt$geneset.names==as.character(gs))
			if(is.null(group.var)){
				if(is.data.frame(expr) | is.matrix(expr)){
					select_probe <- intersect(rownames(expr), unique(gmt$genesets[[interest]]))
					if(!is.numeric(expr)){
						data_sel <- as.matrix(apply(expr[select_probe,], 2, as.numeric))
					}else{
						data_sel <- as.matrix(expr[select_probe,])
					}
				}
				else if(is.list(expr)){
					expr_sel <- expr[[interest]]
					expr_sel <- expr_sel[, , order(as.numeric(dimnames(expr_sel)[[3]]))]
					# data_sel <- matrix(expr_sel, nrow=dim(expr_sel)[1], ncol=dim(expr_sel)[2]*dim(expr_sel)[3])
					data_sel <- acast(melt(expr_sel, varnames=c("Probe_ID", "Subject_ID", "TimePoint")), 
									  formula=Probe_ID ~ TimePoint + Subject_ID)
					# rownames(data_sel) <- dimnames(expr_sel)[[1]]
					select_probe <- dimnames(expr_sel)[[1]]
					TimePoint <- sort(as.numeric(rep(dimnames(expr_sel)[[3]], dim(expr_sel)[2])))
					Subject_ID <- rep(dimnames(expr_sel)[[2]], dim(expr_sel)[3])
				}
				
				data_stand <- t(apply(X=data_sel, MARGIN=1, FUN=scale))
				data_stand[unique(which(is.nan(data_stand), arr.ind=TRUE)[,1]), ] <- 0 
				if(indiv=="genes"){
					data_stand_ByTP <- t(apply(X=data_stand, MARGIN=1, FUN=Fun_byIndex, index=as.factor(TimePoint), fun=aggreg.fun, na.rm=T))
				}
				else if(indiv=="patients"){
					data_tocast<-cbind.data.frame(TimePoint, Subject_ID, "M" = apply(X=data_stand, MARGIN=2, FUN=aggreg.fun, na.rm=T))
					data_stand_ByTP <- as.matrix(acast(data_tocast, formula="Subject_ID~TimePoint", value.var="M"))
				}
				
				if(!is.null(baseline)){
					colbaseline <- which(sort(unique(TimePoint))==baseline)
					if(length(colbaseline)==0){
						stop("the 'baseline' value used is not one of the time points in 'TimePoint'...\n")
					}
					data_stand_ByTP <- data_stand_ByTP-data_stand_ByTP[,colbaseline]
				}
				
			}
			else{
				if(!is.null(baseline)){
					stop("the 'baseline' argument is not NULL while a grouping variable is supplied in 'group.var'...\n")
				}
				if(is.data.frame(expr) | is.matrix(expr)){
					select_probe <- intersect(rownames(expr), unique(gmt$genesets[[interest]]))
					if(!is.numeric(expr)){
						data_sel <- as.matrix(apply(expr[select_probe,], 2, as.numeric))
					}else{
						data_sel <- as.matrix(expr[select_probe,])
					}
				}else if(is.list(expr)){
					expr_sel <- expr[[interest]]
					expr_sel <- expr_sel[, , order(as.numeric(dimnames(expr_sel)[[3]]))]
					data_sel <- matrix(expr_sel, nrow=dim(expr_sel)[1], ncol=dim(expr_sel)[2]*dim(expr_sel)[3])
					rownames(data_sel) <- dimnames(expr_sel)[[1]]
					select_probe <- dimnames(expr_sel)[[1]]
					if(!is.null(Group_ID_paired)){
						Group_ID_paired <- Group_ID_paired[order(TimePoint)] # watch out for the ordering
					}
					TimePoint <- sort(as.numeric(rep(dimnames(expr_sel)[[3]], dim(expr_sel)[2])))
					Subject_ID <- rep(dimnames(expr_sel)[[2]], dim(expr_sel)[3])
				}
				data_stand_ref <- t(apply(X=data_sel[,group.var==ref], MARGIN=1, FUN=scale))
				data_stand_interest <- t(apply(X=data_sel[,group.var==group_of_interest], MARGIN=1, FUN=scale))
				
				if(is.null(Group_ID_paired)){
					data_stand_ByTP_ref <- t(apply(X=data_stand_ref, MARGIN=1, FUN=Fun_byIndex, index=as.factor(TimePoint)[group.var==ref], fun=aggreg.fun))
					data_stand_ByTP_interest <- t(apply(X=data_stand_interest, MARGIN=1, FUN=Fun_byIndex, index=as.factor(TimePoint)[group.var==group_of_interest], fun=aggreg.fun))                         
					data_stand_ByTP <- data_stand_ByTP_interest-data_stand_ByTP_ref
				}else{
					data_stand <- t(apply(X=cbind.data.frame(data_stand_interest, -data_stand_ref), MARGIN=1, FUN=Fun_byIndex, 
										  index=(as.factor(c(TimePoint[group.var==group_of_interest], TimePoint[group.var==ref])):as.factor(c(as.character(Group_ID_paired)[group.var==group_of_interest], as.character(Group_ID_paired)[group.var==ref]))),
										  fun=sum))
					data_stand_ByTP <- t(apply(X=data_stand, MARGIN=1, FUN=Fun_byIndex, index=sort(as.factor(TimePoint[group.var==group_of_interest])), fun=aggreg.fun))
				}
			}
			
			if(sum(apply(data_stand_ByTP, MARGIN=2, FUN=var))<1.e-25){
				nc <- 1
				clust <- rep(1, dim(data_stand)[1])
			} 
			else{
				kmax <- ifelse(dim(data_stand_ByTP)[1]>4, max_trends, dim(data_stand_ByTP)[1]-1)
				if(kmax>=2){
					if(clustering_metric!="sts"){
						cG <- clusGap(x=data_stand_ByTP, FUNcluster=FUNcluster, K.max=kmax, B=B, verbose=FALSE)
						nc <- maxSE(f = cG$Tab[, "gap"], SE.f = cG$Tab[, "SE.sim"], method = methodOptiClust)
						clust <- FUNcluster(data_stand_ByTP, k=nc)$cluster
					}else{
						cG <- clusGap(x=data_stand_ByTP, FUNcluster=FUNcluster, K.max=kmax, B=B, verbose=FALSE, time=as.numeric(colnames(data_stand_ByTP)))
						nc <- maxSE(f = cG$Tab[, "gap"], SE.f = cG$Tab[, "SE.sim"], method = methodOptiClust)
						clust <- FUNcluster(data_stand_ByTP, k=nc, time=as.numeric(colnames(data_stand_ByTP)))$cluster
					}
				}else{
					nc <- 1
					clust <- rep(1, dim(data_stand)[1])
				}
			}
			
			medoids <- as.data.frame(t(apply(X=data_stand_ByTP, MARGIN=2, FUN=Fun_byIndex, index=clust, fun=trend.fun)))
			if(dim(medoids)[1]==1){
				medoids <- cbind.data.frame("TimePoint"= colnames(medoids), "1"=t(medoids))
			}else{
				medoids <- cbind.data.frame("TimePoint"= rownames(medoids), medoids)
			}
			
			NbClust[gs] <- nc
			ClustsMeds[[gs]] <- medoids
			GenesPartition[[gs]] <- clust
			if(indiv=="genes"){
				names(GenesPartition[[gs]]) <- select_probe
			}
			else if(indiv=="patients"){
				names(GenesPartition[[gs]]) <- rownames(data_stand_ByTP)
			}
			if(verbose){
				cat(paste(which(GeneSetsList==gs), "/", length(GeneSetsList), " gene sets clustered\n", sep=""))
			}
		}
		res <- list("NbClust"=NbClust, "ClustMeds"=ClustsMeds, "GenesPartition"=GenesPartition, "MaxNbClust"=max_trends)
		class(res) <- "ClusteredTrends"
		return(res)
	}



#'@rdname clustTrend
#'
#'@param x an object of class '\code{ClusteredTrends}'.
#'
#'@param \dots further arguments passed to or from other methods.
#'
#'@method print ClusteredTrends
#'
#'@export
#'
#'
print.ClusteredTrends <- function(x, ...){
	maxK <- x$MaxNbClust
	f <- factor(x$NbClust, levels=c(1:maxK))
	levels(f)[-1] <- paste(levels(f)[-1], "trends")
	levels(f)[1] <- paste(levels(f)[1], "trend")
	cat("\t\t\tA ClusteredTrends object")
	cat("\n")
	cat("\n")
	cat("Distribution of the number of trends per gene sets:")
	cat("\n\t")
	for(i in 1:maxK){
		cat(levels(f)[i])
		cat(": ")
		if(i==1){cat(" ")}
		cat(summary(f)[i])
		cat("\n\t")
	}
	cat("Total number of trends:", sum(x$NbClust), "(out of", length(x$NbClust), "significant gene sets)\n") 
	cat("\n")
	
	cat("Maximal number of clusters tested:", maxK, "\n")
	
	cat("\n")
	cat("Mean number of trends by significant gene set:", formatC(mean(x$NbClust), digits=3), "\n") 
}



#'@rdname clustTrend
#'
#'@method plot ClusteredTrends
#'
#'@export
#'
plot.ClusteredTrends <- function(x, ...){
	maxK <- x$MaxNbClust
	f <- factor(x$NbClust, levels=c(1:maxK))
	graphics::barplot(height=summary(f),
					  xlab="Number of distinct trends", ylab= "Number of gene sets",
					  col=grDevices::rainbow(maxK),
					  main=paste(formatC(mean(x$NbClust), digits=3), "trends by significant gene set (on average)"),
					  ylim=c(0, sum(summary(f)))  				
	)      
}
