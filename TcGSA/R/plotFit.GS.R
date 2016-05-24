#'Plotting function for exploring the fittness of the mixed modeling used in TcGSA
#'
#'This function plots graphs informing on the fit of the mixed modeling
#'of the gene expression performed in TcGSA, for 1 or several gene sets.
#'
#'@param x a \bold{tcgsa} object for \code{clustTrend}, or a
#'\bold{ClusteredTrends} object for \code{print.ClusteredTrends} and
#'\code{plot.ClusteredTrends}.
#'
#'@param expr 
#'a matrix or dataframe of gene expression.  Its dimension are
#'\eqn{n}x\eqn{p}, with the \eqn{p} samples in column and the \eqn{n} genes in
#'row.
#'
#'@param design
#'a matrix or dataframe containing the experimental variables that used in the model,
#'namely \code{subject_name}, \code{time_name}, and \code{covariates_fixed} 
#'and \code{time_covariates} if applicable.  Its dimension are \eqn{p}x\eqn{m} 
#'and its row are is in the same order as the columns of \code{expr}.
#'
#'@param subject_name
#'the name of the factor variable from \code{design} that contains the information on 
#'the repetition units used in the mixed model, such as the patient identifiers for instance.  
#'Default is \code{'Patient_ID'}.  See Details.
#'
#'@param time_name
#'the name of a numeric variable from \code{design} that contains 
#'the information on the time replicates (the time points at which gene 
#'expression was measured).  Default is \code{'TimePoint'}.  See Details.
#'
#'
#'@param colnames_ID
#'the name of the variable from \code{design} that contains the columnames of the \code{expr}
#'expression data matrix.  See Details.
#'
#'@param plot_type
#'a character string indicating the type of plot to be drawn.  The options are 
#'\code{'Fit'}, \code{'Residuals Obs'}, \code{'Residuals Est'} or \code{'Histogram Obs'}.
#'
#'@param GeneSetsList
#'a character string containing the names of the gene set whose fit is being checked. 
#'If several gene sets are being checked, can be a character list or vector of the 
#'names of those gene sets. 
#'
#'@param color
#'a character string indicating which color scale should be used. One of the 3 : 
#'\code{'genes'}, \code{'time'}, \code{'subjects'}, otherwise, no coloring is used.
#'
#'@param marginal_hist
#'a logical flag indicating wether marginal histograms should be drawn. 
#'Only used for \code{'Fit'} plot type. Default is \code{'TRUE'}
#'
#'@param gg.add 
#'A list of instructions to add to the ggplot2 instruction.  See \link{+.gg}.  Default is \code{list(theme())}, which adds nothing
#'to the plot.
#'
#'@author Boris P. Hejblum
#'
#'@seealso \code{\link{plot1GS}}, \code{\link{plotSelect.GS}}
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
#'@importFrom grDevices rgb
#'
#'@importFrom stats lm
#'
#'@export plotFit.GS
#'
#'@examples
#'data(data_simu_TcGSA)
#'
#'tcgsa_sim_1grp <- TcGSA.LR(expr=expr_1grp, gmt=gmt_sim, design=design, 
#'                           subject_name="Patient_ID", time_name="TimePoint",
#'                           time_func="linear", crossedRandom=FALSE)
#'plotFit.GS(x=tcgsa_sim_1grp, expr=expr_1grp, design=design,
#'					 subject_name="Patient_ID", time_name="TimePoint",
#'					 colnames_ID="Sample_name", 
#'					 plot_type="Residuals Obs", 
#'					 GeneSetsList=c("Gene set 1", "Gene set 2", "Gene set 3",
#'					                "Gene set 4", "Gene set 5"),
#'					 color="genes", gg.add=list(guides(color=FALSE))
#')
#'
#'\dontrun{
#'plotFit.GS(x=tcgsa_sim_1grp, expr=expr_1grp, design=design,
#'           subject_name="Patient_ID", time_name="TimePoint",
#'           colnames_ID="Sample_name", 
#'           plot_type="Histogram Obs", 
#'           GeneSetsList=c("Gene set 1", "Gene set 5"),
#'           color="genes", gg.add=list(guides(fill=FALSE))
#'           )
#'           
#'plotFit.GS(x=tcgsa_sim_1grp, expr=expr_1grp, design=design,
#'           subject_name="Patient_ID", time_name="TimePoint",
#'           colnames_ID="Sample_name", 
#'           plot_type="Histogram Obs", 
#'           GeneSetsList=c("Gene set 1", "Gene set 2", "Gene set 3",
#'			                "Gene set 4", "Gene set 5"),
#'           color="genes")
#'}
#'

plotFit.GS <- function(x, expr, design, subject_name = "Patient_ID", time_name = "TimePoint",
					   colnames_ID, 
					   plot_type=c("Fit", "Residuals Obs", "Residuals Est", "Histogram Obs"),
					   GeneSetsList,
					   color=c("genes", "time", "subjects"),
					   marginal_hist=TRUE,
					   gg.add=list(theme())){
	
	gmt <- x[["GeneSets_gmt"]]
	
	Xfinal <- NULL
	Yfinal <- NULL
	for(gs in GeneSetsList){
		
		interest <- which(gmt$geneset.names==as.character(gs))
		
		select_probe <- intersect(rownames(expr), unique(gmt$genesets[[interest]]))
		if(!is.numeric(expr)){
			xx <- as.matrix(apply(expr[select_probe,], 2, as.numeric))
		}else{
			xx <- as.matrix(expr[select_probe,])
		}
		rownames(xx) <- select_probe
		TimePoint <- design$TimePoint
		Subject_ID <- design$Patient_ID
		xxx <- reshape2::melt(xx, varnames=c("Probe_ID", colnames_ID))
		xxxx <- merge(xxx, design, by=colnames_ID)[, c("Probe_ID", subject_name, time_name, "value")]
		X <- xxxx[order(xxxx$Probe_ID, xxxx[, subject_name], xxxx[, time_name]),]
		X$GS <- gs
		Xfinal <- rbind.data.frame(Xfinal, X)
		
		y <- x[["Estimations"]][[interest]]
		yy <- reshape2::melt(y, varnames=c("Probe_ID", subject_name, time_name))
		Y <- yy[order(yy$Probe_ID, yy[ ,subject_name], yy[, time_name]),]
		Y$GS <- gs
		Yfinal <- rbind.data.frame(Yfinal, Y)
	}
	
	data2plot <- merge(Xfinal,Yfinal, by=c("Probe_ID", subject_name, time_name, "GS"))
	colnames(data2plot) <- c("Probe_ID", "Subject", "Time", "GS", "Observed", "Estimated")
	data2plot[, "Time"] <- as.factor(data2plot[, "Time"])
	data2plot[, "Probe_ID"] <- as.factor(data2plot[, "Probe_ID"])
	data2plot$Residuals <- data2plot$Observed - data2plot$Estimated
	
	
	if(plot_type=="Histogram Obs"){
		p <- (ggplot(aes_string(x="Observed"), data=data2plot)
			  + labs(x="Observed expression", y="Count")
			  + theme_bw()
			  
		)
		if(color=="time"){
			p <- (p + geom_histogram(aes_string(fill="Time"))
				  + scale_fill_discrete(name="Time point")
				  + ggtitle(paste(gs, ": ", gmt$geneset.descriptions[interest],"\n Histogram of observed expression", sep=""))
			)
		}else if(color=="genes"){
			p <- (p + geom_histogram(aes_string(fill="Probe_ID"))			
				  + scale_fill_discrete(name="Genes")
				  + ggtitle(paste(gs, ": ", gmt$geneset.descriptions[interest],"\n Histogram of observed expression colored by genes", sep=""))
			)
		}else if(color=="subjects"){
			p <- (p + geom_histogram(aes_string(fill="Subject"))			
				  + scale_fill_discrete(name="Subjects")
				  + ggtitle(paste(gs, ": ", gmt$geneset.descriptions[interest],"\n Histogram of observed expression", sep=""))
			)
		}
		else{
			warning("'color' argument must be either 'genes', 'time', 
					or 'subject' when 'Histogram Obs' are plotted")
			p <- p + geom_histogram()
		}
	}else{
		if(plot_type=="Fit"){
			plotmin <- min(c(data2plot$Estimated, data2plot$Observed))
			plotmax <- max(c(data2plot$Estimated, data2plot$Observed))
			p <- (ggplot(aes_string(x="Observed", y="Estimated"), data="data2plot") 
				  + geom_smooth(method=lm)
				  + labs(title=paste(gs, ": ", gmt$geneset.descriptions[interest],"\n Fit", sep=""),
				  	   x="Observed expression", y="Fitted expression")
				  + xlim(plotmin, plotmax)
				  + ylim(plotmin, plotmax)
				  + theme_bw()	
			)
			if(color=="genes"){
				p <- p + ggtitle(paste(gs, ": ", gmt$geneset.descriptions[interest],"\n Fit colored by genes", sep=""))
			}
		}else if(plot_type=="Residuals Obs"){
			p <- (ggplot(aes_string(x="Observed", y="Residuals"), data=data2plot)
				  + labs(title=paste(gs, ": ", gmt$geneset.descriptions[interest],"\n Residuals", sep=""),
				  	   x="Observed expression", y="Residuals")
				  #+ xlim(min(data2plot$Observed),max(data2plot$Observed))
				  #+ ylim(-max(abs(data2plot$Residuals)), max(abs(data2plot$Residuals)))
				  + theme_bw()
			)
			if(color=="genes"){
				p <- p + ggtitle(paste(gs, ": ", gmt$geneset.descriptions[interest],"\n Residuals colored by genes", sep=""))
			}
		}else if(plot_type=="Residuals Est"){
			p <- (ggplot(aes_string(x="Estimated", y="Residuals"), data=data2plot)
				  + labs(title=paste(gs, ": ", gmt$geneset.descriptions[interest],"\n Residuals", sep=""),
				  	   x="Estimated expression", y="Residuals")
				  #+ xlim(min(data2plot$Estimated),max(data2plot$Estimated))
				  #+ ylim(-max(abs(data2plot$Residuals)), max(abs(data2plot$Residuals)))
				  + theme_bw()	
			)
			if(color=="genes"){
				p <- p + ggtitle(paste(gs, ": ", gmt$geneset.descriptions[interest],"\n Residuals colored by genes", sep=""))
			}
		}
		if(color=="time"){
			p <- (p + geom_point(aes_string(color="Time"))
				  + scale_color_discrete(name="Time point")
			)
		}else if(color=="genes"){
			p <- (p + geom_point(aes_string(color="Probe_ID"))
				  + scale_color_discrete(name="Genes")
			)
		}else if(color=="subjects"){
			p <- (p + geom_point(aes_string(color="Subject"))
				  + scale_color_discrete(name="Subjects")
			)
		}else{
			warning("wrong 'color' argument")
			p <- p + geom_point()
		}
	}
	
	if(length(GeneSetsList)>1){
		p <- (p + facet_wrap( ~GS, ncol=floor(sqrt(length(unique(data2plot$GS)))),
							  scales="free")
			  + ggtitle(plot_type)
		)
	}
	
	if(plot_type=="Fit"){
		p <- p + geom_abline(intercept=0, slope=1, color="black")
		if(marginal_hist){
			p <- p + geom_rug(, col=grDevices::rgb(0,0,0,alpha=.3))
		}
	}else if(plot_type=="Residuals Obs" | plot_type=="Residuals Est" ){
		p <- (p + geom_abline(intercept=0, slope=0, color="red")
			  + geom_smooth(method=lm))
		
	}
	
	for(a in gg.add){
		p <- p + a
	}
	
	return(p)
	
}



