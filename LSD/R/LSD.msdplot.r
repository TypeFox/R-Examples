

### singlemsdplot ###


#' @export
#' @name singlemsdplot
#' @aliases LSD.singlemsdplot
#' @title Visualize two-dimensional data clusters (add to an existing plot)
#' @description Depict a numeric matrix or list utilizing the underlying mean and standard deviation estimates of one dimension in a color encoded fashion (add to an existing plot).
#' @param input data as matrix or list.
#' @param col a character vector of R build-in colors.
#' @param alpha alpha value: a two-digit integer between 01 and 99 for color opacity, i.e. appearance of partial or full transparency (usage omitted by default).
#' @param bars logical: if \code{TRUE} (by default), error bars are added at each position.
#' @param length a numeric value scaling the width of the bars.
#' @param at a integer vector containing the x-positions corresponding to the columns of 'input'.
#' @author Bjoern Schwalb
#' @seealso \code{\link{comparisonplot}}, \code{\link{demotour}}, \code{\link{disco}}, \code{\link{colorpalette}}
#' @examples samples = 100
#' probes = 200
#' clus = matrix(rnorm(probes*samples,sd=1),ncol=probes)
#' 
#' clus = rbind(
#' 	t(t(clus)+sin(1:probes/10))+1:nrow(clus)/samples,
#' 	t(t(clus)+sin(pi/2+1:probes/10))+1:nrow(clus)/samples)
#' 
#' emptyplot(xlim = c(1,ncol(clus)),ylim = range(clus))
#' singlemsdplot(clus)
#' axis(1)
#' axis(2)
#' box()
#' @keywords mean, sd


singlemsdplot = function(input,col = "darkgreen",alpha = 50,bars = TRUE,length = 0.25,at = NULL)
{
	if (!is.matrix(input) & !is.list(input)) stop("First argument must be a matrix or a list!")
	if (is.null(at)) if (is.matrix(input)){at=1:ncol(input)} else if (is.list(input)){at=1:length(input)}
	probes = length(at)
	if (is.matrix(input)){
		means = apply(input,2,mean,na.rm=TRUE)
		sds = apply(input,2,sd,na.rm=TRUE)
	} else if (is.list(input)){
		means = lapply(input,mean,na.rm=TRUE)
		means = sapply(means,c)
		sds = lapply(input,sd,na.rm=TRUE)
		sds = sapply(sds,c)
	}
	sds[is.na(sds)] = 0
	qline = rbind(means+sds,means,means-sds)
	for (j in 1:2){polygon(at[c(1:probes,probes:1)],c(qline[j,],qline[j+1,probes:1]),col = convertcolor(col,alpha),lty=0)}
	lines(at,qline[2,],col=col,lwd=2)
	if (bars){options(warn = -1);arrows(at,means-sds,at,means+sds,angle=90,code=3,length=length,col=col);options(warn = 0)}                               
}


### aliases ###


LSD.singlemsdplot = singlemsdplot


### msdplot ###


#' @export
#' @name msdplot
#' @aliases LSD.msdplot
#' @title Visualize two-dimensional data clusters
#' @description Depict a numeric matrix or list utilizing the underlying mean and standard deviation estimates of one dimension in a color encoded fashion.
#' @param input matrix or list with numerical entries, quantiles of cols will define lines.
#' @param label a character vector assigning rows/elements of 'input' to clusters (if specified, multiple clusters can be depicted in different colors and/or subsequent plots).
#' @param at a integer vector containing the x-positions corresponding to the columns of 'input'.
#' @param xlim x limits, standard graphics parameter.
#' @param ylim y limits, standard graphics parameter.
#' @param xlab x labels, standard graphics parameter.
#' @param ylab y labels, standard graphics parameter.
#' @param main title(s) of the plot, standard graphics parameter.
#' @param xaxt a character which specifies the x axis type ("n" suppresses plotting of the axis).
#' @param xlabels a character vector containing labels for the x-axis.
#' @param las las=1: horizontal text, las=2: vertical text (x-axis labels).
#' @param separate if \code{TRUE} (by default), different clusters are depicted in subsequent plots.
#' @param size logical: if \code{TRUE} (by default), the size of each cluster is added to the title of the respective plot.
#' @param col a character vector giving R build-in colors for different clusters.
#' @param bars logical: if \code{TRUE} (by default), error bars are added at each position.
#' @param alpha alpha value: a two-digit integer between 01 and 99 for color opacity, i.e. appearance of partial or full transparency (usage omitted by default).
#' @param ... additional parameters to be passed to points and plot.
#' @author Bjoern Schwalb
#' @seealso \code{\link{comparisonplot}}, \code{\link{demotour}}, \code{\link{disco}}, \code{\link{colorpalette}}
#' @examples at = c(2,4,8,16,32)
#' clus = matrix(rnorm(500,sd=0.5),ncol=5)
#' batch = sample(c(-8,-6,-4,-2),100,replace=TRUE)
#' clus = clus + cbind(0,0.25*batch,0.5*batch,0.75*batch,batch)
#' clus = clus - clus[,1]
#' clus = t(t(clus)*c(0,0.1,0.25,0.5,1))
#' labs = paste("cluster",kmeans(clus,4)$cluster)
#' 
#' colpal = c("darkgreen","darkblue","darkred","black")
#' msdplot(clus,labs,at,separate=FALSE,col=colpal,alpha=25,xlabels=at)
#' 
#' msdplot(clus,labs,at,col=colpal,alpha=50,xlabels=at)
#' @keywords mean, sd


msdplot = function(input,label = NULL,at = NULL,xlim = NULL,ylim = NULL,xlab = "",ylab = "",main = "msdplot",xaxt = "s",xlabels = NULL,las = 1,separate = TRUE,size = TRUE,col = "darkgreen",bars = TRUE,alpha = 50,...)
{
	if (!is.matrix(input) & !is.list(input)) stop("First argument must be a matrix or a list!")
	if (is.null(at)) if (is.matrix(input)){at=1:ncol(input)} else if (is.list(input)){at=1:length(input)}
	probes = length(at)
	if (is.null(xlim)) xlim=range(at)
	maxp = xlim[2]
	minp = xlim[1]
	if (is.null(ylim)) ylim=range(input,na.rm=TRUE)
	if (is.null(xlabels)) xlabels = 1:length(at)
	# one cluster in one plot
	if (is.null(label)){
		plot(xlim,ylim,type="n",xaxt="n",xlab=xlab,ylab=ylab,main="",...)
		mtext(paste(main),3,2,cex=1.25)
		if (size){mtext(paste("( #",nrow(input)," )"),3,1,col="darkgrey",cex=0.75)}
		if (xaxt!="n") axis(side=1,las=las,at=at,labels=xlabels)
		singlemsdplot(input=input,at=at,alpha=alpha,col=col,bars=bars,length=0.02*range(at))
	}
	# several clusters in one plot or several plots
	if (!is.null(label)){clusternames = unique(label)
		nrclusters = length(clusternames)
		if (!is.matrix(input)) stop("First argument must be a matrix for multiple clusters!")
		clustersets = split(1:nrow(input),factor(label))	
		if (separate==FALSE){
			plot(xlim,ylim,type="n",xaxt="n",xlab=xlab,ylab=ylab,main="",...)
			mtext(paste(main),3,2,cex=1.25)
			if (size){mtext(paste("( #",nrow(input)," )"),3,1,col="darkgrey",cex=0.75)}
			if (xaxt!="n") axis(side=1,las=las,at=at,labels=xlabels)}
		if (separate==TRUE) par(mfrow=windowxy(nrclusters))
		for (j in seq(clusternames)){
			if (separate==TRUE){
				if (length(main) == length(clustersets[[j]])) clustermain = main[j]	else clustermain = paste(main,clusternames[j])
				plot(xlim,ylim,type="n",xaxt="n",xlab=xlab,ylab=ylab,main="",...)
				mtext(paste(clustermain),3,2,cex=1.25)
				if (size){mtext(paste("( #",length(clustersets[[j]])," )"),3,1,col="darkgrey",cex=0.75)}
				if (xaxt!="n") axis(side=1,las=las,at=at,labels=xlabels)
			}
			singlemsdplot(input=input[clustersets[[j]],,drop=FALSE],at=at,alpha=alpha,col=col[j],bars=bars,length=0.02*range(at))
		}
	}
}


### aliases ###


LSD.msdplot = msdplot



