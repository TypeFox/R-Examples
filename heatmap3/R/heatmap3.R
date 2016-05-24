##' heatmap3
##' 
##' The function heatmap3 is completely compatible with the original R function heatmap, and provides more new features.
##' 
##' 
##' 
##' @inheritParams stats::heatmap
##' @param useRaster logical; if TRUE a bitmap raster is used to plot the image instead of polygons. The grid must be regular in that case, otherwise an error is raised.
##' @param file pdf file name, only works when topN was used.
##' @param topN vector a list of numbers. topN genes will be used to generate the heatmaps.
##' @param filterFun function used to filter genes, such as sd, mean, sum. It will be used in a apply function to caculate for each row.
##' @param legendfun function used to generate legend in top left of the figure. If not specified, the color bar will be plotted. The users can use any plot functions to generate their own legend. Or a function \code{\link{showLegend}} is also provided as a example.
##' @param ColSideFun function used to generate annotation and labeling figure in column side. The users can use any plot functions to generate their own figure. And a function \code{\link{showAnn}} is also provided as a example.
##' @param ColSideAnn data frame with continuous and factor variables as annotation information. This parameter will be sorted by coloum dendrogram and then passed to ColSideFun.
##' @param ColSideWidth numeric the height of column side area, which can be used by ColSideFun function.
##' @param ColSideCut numeric the value to be used in cutting coloum dendrogram. The dendrogram and annotation will be divided into different parts and labeled respectively.
##' @param colorCell A data.frame with 3 columns, indicating which cells will be colored by specific colors. The first column is row index, second column is column index, and the third column is color.
##' @param highlightCell A data.frame with 3 or 4 columns, indicating which cells will be highlighted by rectangles with specific colors. The first column is row index, second column is column index, the third column is color for rectangle border, and the optional forth column is width for rectangle border. 
##' @param method the agglomeration method to be used by \code{\link{hclust}} function. This should be (an unambiguous abbreviation of) one of "ward", "single", "complete", "average", "mcquitty", "median" or "centroid".
##' @param balanceColor logical indicating if the colors need to be balanced so that the median color will represent the 0 value. The default value is F.
##' @param ColAxisColors integer indicating which coloum of ColSideColors will be used as colors for labels in coloum axis. The default value is 0, which means all coloum labels will be in black color.
##' @param RowAxisColors integer indicating which coloum of RowSideColors will be used as colors for labels in row axis. The default value is 0, which means all row labels will be in black color.
##' @param showColDendro logical indicating if the coloum dendrogram should be plotted (when Colv isn't NA).
##' @param showRowDendro logical indicating if the row dendrogram should be plotted (when Rowv isn't NA).
##' @param RowSideLabs label for RowSideColors
##' @param ColSideLabs label for ColSideColors
##' @param col specifying the colors, used in \code{\link{image}} function.
##' @param cexRow,cexCol positive numbers, used as cex.axis in for the row or column axis labeling. The defaults currently only use number of rows or columns, respectively.
##' @param labRow,labCol character vectors with row and column labels to use; these default to rownames(x) or colnames(x), respectively.
##' @param lasRow,lasCol the style of row or column axis labels.
##' @param main,xlab,ylab main, x- and y-axis titles; defaults to none.
##' @param ... additional arguments passed on to \code{\link{image}}.
##' @importFrom fastcluster hclust
##' @export
##' @return The same return value as \code{\link{hclust}} function.
##' @examples #gererate data
##' set.seed(123456789)
##' rnormData<-matrix(rnorm(1000), 40, 25)
##' rnormData[1:15, seq(6, 25, 2)] = rnormData[1:15, seq(6, 25, 2)] + 2
##' rnormData[16:40, seq(7, 25, 2)] = rnormData[16:40, seq(7, 25, 2)] + 4
##' colnames(rnormData)<-c(paste("Control", 1:5, sep = ""), paste(c("TrtA", "TrtB"),
##' rep(1:10,each=2), sep = ""))
##' rownames(rnormData)<-paste("Probe", 1:40, sep = "")
##' ColSideColors<-cbind(Group1=c(rep("steelblue2",5), rep(c("brown1", "mediumpurple2"),10)),
##'     Group2=sample(c("steelblue2","brown1", "mediumpurple2"),25,replace=TRUE))
##' colorCell<-data.frame(row=c(1,3,5),col=c(2,4,6),color=c("green4","black","orange2"),
##'     stringsAsFactors=FALSE)
##' highlightCell<-data.frame(row=c(2,4,6),col=c(1,3,5),color=c("black","green4","orange2"),
##'     lwd=1:3,stringsAsFactors=FALSE)
##' #A simple example
##' heatmap3(rnormData,ColSideColors=ColSideColors,showRowDendro=FALSE,colorCell=colorCell,
##'     highlightCell=highlightCell)
##' #A more detail example
##' ColSideAnn<-data.frame(Information=rnorm(25),Group=c(rep("Control",5), rep(c("TrtA", "TrtB"),10)))
##' row.names(ColSideAnn)<-colnames(rnormData)
##' RowSideColors<-colorRampPalette(c("chartreuse4", "white", "firebrick"))(40)
##' result<-heatmap3(rnormData,ColSideCut=1.2,ColSideAnn=ColSideAnn,ColSideFun=function(x) 
##' showAnn(x),ColSideWidth=0.8,RowSideColors=RowSideColors,col=colorRampPalette(c("green",
##' "black", "red"))(1024),RowAxisColors=1,legendfun=function() showLegend(legend=c("Low",
##' "High"),col=c("chartreuse4","firebrick")),verbose=TRUE)
##' #annotations distribution in different clusters and the result of statistic tests
##' result$cutTable
heatmap3<-function (x, Rowv = NULL, Colv = if (symm) "Rowv" else NULL, 
		distfun = function(x) as.dist(1 - cor(t(x),use="pa")),balanceColor=F, ColSideLabs,RowSideLabs,showColDendro=T,showRowDendro=T,col=colorRampPalette(c("navy", "white", "firebrick3"))(1024),legendfun,method="complete",ColAxisColors=0,RowAxisColors=0, hclustfun = hclust, reorderfun = function(d, 
				w) reorder(d, w), add.expr,symm = FALSE, revC = identical(Colv, 
				"Rowv"), scale = c("row", "column", "none"), na.rm = TRUE, 
		ColSideFun,ColSideAnn,ColSideWidth=0.4,ColSideCut,colorCell,highlightCell,
		file="heatmap3.pdf",topN=NA,filterFun=sd,
		margins = c(5, 5), ColSideColors, RowSideColors, cexRow = 0.2 + 
				1/log10(nrow(x)), cexCol = 0.2 + 1/log10(ncol(x)),lasRow=2,lasCol=2, labRow = NULL, 
		labCol = NULL, main = NULL, xlab = NULL, ylab = NULL, keep.dendro = FALSE, 
		verbose = getOption("verbose"),useRaster=if (ncol(x)*nrow(x)>=50000) TRUE else FALSE ,...) 
{
	#loop fot different topN
	if (!all(is.na(topN))) {
		temp<-apply(x,1,filterFun)
		pdf(file)
		for (n in topN) {
			xSub<-x[rev(order(temp))[1:n],,drop=F]
			if (!missing(RowSideColors)) {
				RowSideColorsBak<-RowSideColors
				RowSideColors<-RowSideColors[rev(order(temp))[1:n],,drop=F]
			}
			result[[paste0(n)]]<-heatmap3(xSub,Rowv=Rowv,Colv=Colv,distfun=distfun,balanceColor=balanceColor,ColSideLabs=ColSideLabs,RowSideLabs=RowSideLabs,showColDendro=showColDendro,showRowDendro=showRowDendro,col=col,legendfun=legendfun,method="complete",ColAxisColors=0,RowAxisColors=0, hclustfun = hclust, reorderfun=reorderfun,
					add.expr=add.expr,symm = symm, revC=revC,scale=scale,na.rm=na.rm,
					ColSideFun=ColSideFun,ColSideAnn=ColSideAnn,ColSideWidth=ColSideWidth,ColSideCut=ColSideCut,
					margins = margins,ColSideColors=ColSideColors,RowSideColors=RowSideColors,cexRow = cexRow,
					cexCol = cexCol, labRow = labRow, 
					labCol = labCol, main = paste0("top ",n), xlab = xlab, ylab = ylab, keep.dendro = keep.dendro, 
					verbose = verbose,...		
			)
			if (!missing(RowSideColors)) {
				RowSideColors<-RowSideColorsBak
			}
		}
		temp<-dev.off()
		cat(paste0("The heatmaps were generated at ",file,"\n"))
		return(invisible(result))
	}
	
	scale <- if (symm && missing(scale)) 
				"none"
			else match.arg(scale)
	if (is.data.frame(x)) {x<-as.matrix(x)}
	if (!missing(ColSideColors)) {
		if (is.vector(ColSideColors)) {
			ColSideColors<-cbind(ColSideColors)
		}
	}
	if (!missing(RowSideColors)) {
		if (is.vector(RowSideColors)) {
			RowSideColors<-cbind(RowSideColors)
		}
	}
	if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
		stop("'x' must be a numeric matrix")
	nr <- di[1L]
	nc <- di[2L]
	if (nr <= 1 || nc <= 1) 
		stop("'x' must have at least 2 rows and 2 columns")
	if (!is.numeric(margins) || length(margins) != 2L) 
		stop("'margins' must be a numeric vector of length 2")
	doRdend <- !identical(Rowv, NA)
	doCdend <- !identical(Colv, NA)
	if (!doRdend && identical(Colv, "Rowv")) 
		doCdend <- FALSE
	if (is.null(Rowv)) 
		Rowv <- rowMeans(x, na.rm = na.rm)
	if (is.null(Colv)) 
		Colv <- colMeans(x, na.rm = na.rm)
	if (doRdend) {
		if (inherits(Rowv, "dendrogram")) 
			ddr <- Rowv
		else {
			hcr <- hclustfun(distfun(x),method=method)
			ddr <- as.dendrogram(hcr)
			if (!is.logical(Rowv) || Rowv) 
				ddr <- reorderfun(ddr, Rowv)
		}
		if (nr != length(rowInd <- order.dendrogram(ddr))) 
			stop("row dendrogram ordering gave index of wrong length")
	}
	else rowInd <- 1L:nr
	if (doCdend) {
		if (inherits(Colv, "dendrogram")) 
			ddc <- Colv
		else if (identical(Colv, "Rowv")) {
			if (nr != nc) 
				stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
			ddc <- ddr
		}
		else {
			hcc <- hclustfun(distfun(if (symm) 
										x
									else t(x)),method=method)
			ddc <- as.dendrogram(hcc)
			if (!is.logical(Colv) || Colv) 
				ddc <- reorderfun(ddc, Colv)
		}
		if (nc != length(colInd <- order.dendrogram(ddc))) 
			stop("column dendrogram ordering gave index of wrong length")
	}
	else colInd <- 1L:nc
	x <- x[rowInd, colInd]
	labRow <- if (is.null(labRow)) 
				if (is.null(rownames(x))) 
					(1L:nr)[rowInd]
				else rownames(x)
			else labRow[rowInd]
	labCol <- if (is.null(labCol)) 
				if (is.null(colnames(x))) 
					(1L:nc)[colInd]
				else colnames(x)
			else labCol[colInd]
	if (scale == "row") {
		x <- sweep(x, 1L, rowMeans(x, na.rm = na.rm), check.margin = FALSE)
		sx <- apply(x, 1L, sd, na.rm = na.rm)
		x <- sweep(x, 1L, sx, "/", check.margin = FALSE)
	}
	else if (scale == "column") {
		x <- sweep(x, 2L, colMeans(x, na.rm = na.rm), check.margin = FALSE)
		sx <- apply(x, 2L, sd, na.rm = na.rm)
		x <- sweep(x, 2L, sx, "/", check.margin = FALSE)
	}
	lmat <- rbind(c(NA, 3), 2:1)
	lwid <- c(1, 4)
	lhei <- c( 1 + if (!is.null(main)) 0.2 else 0,4)
	if (!missing(ColSideFun)) {
		lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
		lhei <- c(lhei[1L], ColSideWidth, lhei[2L])
	} else if (!missing(ColSideColors)){ #ColSideColors
#		if (!is.character(ColSideColors) || nrow(ColSideColors) != 
		if (!is.character(ColSideColors) & nrow(ColSideColors) != 
				nc) 
			stop("'ColSideColors' must be a character vector or matrix of length ncol(x)")
		lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
		lhei <- c(lhei[1L], 0.2*round(ncol(ColSideColors)/2+0.1), lhei[2L])
	}
	if (!missing(RowSideColors)) {
		if (!is.character(RowSideColors) || nrow(RowSideColors) != 
				nr) 
			stop("'RowSideColors' must be a character vector or matrix of length nrow(x)")
		lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 
						1), lmat[, 2] + 1)
		lwid <- c(lwid[1L], 0.2*round(ncol(RowSideColors)/2+0.1), lwid[2L])
	}
	lmat<-lmat+1
	lmat[is.na(lmat)] <- 0
	lmat[1,1]<-1
#	if (verbose) {
#		cat("layout: widths = ", lwid, ", heights = ", lhei, 
#				"; lmat=\n")
#		print(lmat)
#	}
	dev.hold()
	on.exit(dev.flush())
	op <- par(no.readonly = TRUE)
	on.exit(par(op), add = TRUE)
	
	#balanceColor
	if (balanceColor) {
		if (abs(max(x,na.rm=T))>=abs(min(x,na.rm=T))) {
			cut.off<-round(quantile(1:length(col),probs=1-(abs(max(x,na.rm=T))+abs(min(x,na.rm=T)))/(2*abs(max(x,na.rm=T)))))
			col<-col[cut.off:length(col)]
		} else {
			cut.off<-round(quantile(1:length(col),probs=(abs(max(x,na.rm=T))+abs(min(x,na.rm=T)))/(2*abs(min(x,na.rm=T)))))
			col<-col[1:cut.off]
		}
	}
	
	layout(lmat, widths = lwid, heights = lhei, respect = TRUE)
	if (!missing(legendfun)) {
		par(mar = c(0, 0, 0, 0))
		legendfun()
	} else {
		par(mar = c(5, 1, 1, 0))
		dummy.x <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), 
				length = length(col))
		dummy.z <- matrix(dummy.x, ncol = 1)
		image(x = dummy.x, y = 1, z = dummy.z, yaxt = "n", 
				col = col,cex.axis=cexCol,xlab="")
	}
	if (!missing(RowSideColors)) {
		par(mar = c(margins[1L], 0, 0, 0.5))
		if (revC) {
			rsc = RowSideColors[rev(rowInd),,drop=F]
		} else {
			rsc = RowSideColors[rowInd,,drop=F]
		}
		rsc.colors = matrix()
		rsc.names = names(table(rsc))
		rsc.i = 1
		for (rsc.name in rsc.names) {
			rsc.colors[rsc.i] = rsc.name
			rsc[rsc == rsc.name] = rsc.i
			rsc.i = rsc.i + 1
		}
		rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
		image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
		
		if (missing(RowSideLabs)) {
			if (ncol(RowSideColors)==1 & colnames(RowSideColors)[1]=="") {
				RowSideLabs<-""
			} else {
				RowSideLabs<-colnames(RowSideColors)
			}
		}
		if (dim(rsc)[2]==1) {
			axis(1, 0, RowSideLabs, las = 2, tick = FALSE)
		} else {
			axis(1, 0:(dim(rsc)[2] - 1)/(dim(rsc)[2] - 1), RowSideLabs, 
					las = 2, tick = FALSE)
		}
	}
	if (!missing(ColSideCut)) {
		ColSideCutResult<-cut(ddc,ColSideCut)$lower
		cutTable<-NULL
		if (verbose) {
			cat(paste0("The samples could be cut into ",length(ColSideCutResult)," parts with height ",ColSideCut))
			cat("\n")
			ColSideCutResultSubIndList<-list()
			for (i in 1:length(ColSideCutResult)) {
				ColSideCutResultSubInd<-order.dendrogram(ColSideCutResult[[i]])
				ColSideCutResultSubIndList[[i]]<-ColSideCutResultSubInd
#				cat(paste0("Summary for cluster ",i,":\n"))
#				print(summary(ColSideAnn[ColSideCutResultSubInd,]))
			}
			for (i in 1:ncol(ColSideAnn)) {
				if (is.factor(ColSideAnn[,i])) { #factor
					cutTable[[i]]<-sapply(ColSideCutResultSubIndList,function(x) table(ColSideAnn[x,i]))
					colnames(cutTable[[i]])<-paste0("Cluster ",1:length(ColSideCutResult))
					names(cutTable)[i]<-colnames(ColSideAnn)[i]
					pvalue<-chisq.test(cutTable[[i]])$p.value
					cat(paste0("Differential distribution for ",colnames(ColSideAnn)[i],", p value by chi-squared test: ",round(pvalue,3),"\n"))
					cutTable[[i]]<-rbind(cutTable[[i]],round(cutTable[[i]][1,]/colSums(cutTable[[i]]),2))
					row.names(cutTable[[i]])[nrow(cutTable[[i]])]<-paste0(row.names(cutTable[[i]])[1],"_Percent")
					cutTable[[i]]<-cbind(cutTable[[i]],pValue=c(pvalue,rep(NA,nrow(cutTable[[i]])-1)))
				} else { #continous
					cutTable[[i]]<-sapply(split(ColSideAnn[unlist(ColSideCutResultSubIndList),i],rep(1:length(ColSideCutResultSubIndList),sapply(ColSideCutResultSubIndList,length))),function(x) summary(na.omit(x)))
					colnames(cutTable[[i]])<-paste0("Cluster ",1:length(ColSideCutResult))
					names(cutTable)[i]<-colnames(ColSideAnn)[i]
					temp<-aov(ColSideAnn[unlist(ColSideCutResultSubIndList),i]~as.factor(rep(1:length(ColSideCutResultSubIndList),sapply(ColSideCutResultSubIndList,length))))
					pvalue<-summary(temp)[[1]]$"Pr(>F)"[1]
					cat(paste0("Differential distribution for ",colnames(ColSideAnn)[i],", p value by ANOVA: ",round(pvalue,3),"\n"))
					cutTable[[i]]<-cbind(cutTable[[i]],pValue=c(pvalue,rep(NA,5)))
				}
			}
			
		}
		ColSideCutResultCol<-rainbow(length(ColSideCutResult),alpha=0.2)
		ColNumber<-(ncol(x)-1)
	}
	if (!missing(ColSideFun)) {
		par(mar = c(0.5, 0, 0, margins[2L]))
		ColSideAnn<-ColSideAnn[colInd,,drop=F]
		ColAnnHeight<-ColSideFun(ColSideAnn)
		if (!exists("ColAnnHeight")) {
			ColAnnHeight<-par("usr")[3:4]
		}
		if (!missing(ColSideCut)) {
			rect(c(0-1/ColNumber/2,(0-1/ColNumber/2)+1/ColNumber*cumsum(sapply(ColSideCutResult,function(x) length(unlist(x))))[-length(ColSideCutResult)]),ColAnnHeight[1],c((0-1/ColNumber/2)+1/ColNumber*cumsum(sapply(ColSideCutResult,function(x) length(unlist(x))))), ColAnnHeight[2],col=ColSideCutResultCol)
		}
	} else if (!missing(ColSideColors)) {
		par(mar = c(0.5, 0, 0, margins[2L]))
		csc = ColSideColors[colInd,,drop=F]
		csc.colors = matrix()
		csc.names = names(table(csc))
		csc.i = 1
		for (csc.name in csc.names) {
			csc.colors[csc.i] = csc.name
			csc[csc == csc.name] = csc.i
			csc.i = csc.i + 1
		}
		csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
		image(csc, col = as.vector(csc.colors), axes = FALSE)
		if (missing(ColSideLabs)) {
			if (ncol(ColSideColors)==1 & colnames(ColSideColors)[1]=="") {
				ColSideLabs<-""
			} else {
				ColSideLabs<-colnames(ColSideColors)
			}
		}
		if (dim(csc)[2]==1) {
			axis(4, 0, ColSideLabs,las = 2, tick = FALSE)
		} else {
			axis(4, 0:(dim(csc)[2] - 1)/(dim(csc)[2] - 1), ColSideLabs,las = 2, tick = FALSE)
		}
	}
	par(mar = c(margins[1L], 0, 0, margins[2L]))
	if (!symm || scale != "none") 
		x <- t(x)
	if (revC) {
		iy <- nr:1
		if (doRdend)
			ddr <- rev(ddr)
		x <- x[, iy]
	}
	else iy <- 1L:nr
	image(1L:nc, 1L:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
					c(0, nr), axes = FALSE, xlab = "", ylab = "", col=col,useRaster=useRaster,...)
	if (!missing(colorCell)) {
		colorCell[,1]<-rowInd[colorCell[,1]]
		colorCell[,2]<-colInd[colorCell[,2]]
		rect(colorCell[,2]-0.5,colorCell[,1]-0.5,colorCell[,2]+0.5,colorCell[,1]+0.5,col=as.character(colorCell[,3]),border=NA)
	}
	if (!missing(highlightCell)) {
		if (ncol(highlightCell)==3) {
			highlightCell$lwd<-1
		}
		highlightCell[,1]<-rowInd[highlightCell[,1]]
		highlightCell[,2]<-colInd[highlightCell[,2]]
		rect(highlightCell[,2]-0.5,highlightCell[,1]-0.5,highlightCell[,2]+0.5,highlightCell[,1]+0.5,border=as.character(highlightCell[,3]),lwd=as.integer(highlightCell[,4]))
	}
	if (!missing(ColSideColors) & ColAxisColors!=0) {
		mtext(1, at=1L:nc, text = labCol, las = lasCol, line = 0.5,cex = cexCol,col=ColSideColors[colInd,ColAxisColors])
	} else {
		axis(1, 1L:nc, labels = labCol, las = lasCol, line = -0.5, tick = 0,cex.axis = cexCol)
	}
	
	if (!is.null(xlab)) 
		mtext(xlab, side = 1, line = margins[1L] - 1.25)
	if (!missing(RowSideColors) & RowAxisColors!=0) {
		mtext(4, at=iy, text = labRow, las = lasRow, line = 0.5,cex = cexRow,col=RowSideColors[rowInd,RowAxisColors])
	} else {
		axis(4, iy, labels = labRow, las = lasRow, line = -0.5, tick = 0, cex.axis = cexRow)
	}
	if (!is.null(ylab)) 
		mtext(ylab, side = 4, line = margins[2L] - 1.25)
	if (!missing(add.expr)) 
		eval(substitute(add.expr))
	par(mar = c(margins[1L], 0, 0, 0))
	if (doRdend & showRowDendro) 
		plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
	else frame()
	par(mar = c(0, 0, if (!is.null(main)) 1 else 0, margins[2L]))
	if (doCdend & showColDendro) {
		plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
		if (!missing(ColSideCut)) {
			rect(c(0.5,0.5+cumsum(sapply(ColSideCutResult,function(x) length(unlist(x))))[-length(ColSideCutResult)]), 0, cumsum(sapply(ColSideCutResult,function(x) length(unlist(x))))+0.5, ColSideCut,col=ColSideCutResultCol)
		}
	}	else if (!is.null(main)) 
		frame()
	if (!is.null(main)) {
		par(xpd = NA)
		title(main, cex.main = 1.5 * op[["cex.main"]])
	}
	invisible(list(rowInd = rowInd, colInd = colInd, Rowv = if (keep.dendro && 
							doRdend) ddr, Colv = if (keep.dendro && doCdend) ddc, cutTable = if (!missing(ColSideAnn) && !missing(ColSideCut)) cutTable))
}

##' showLegend
##' 
##' The function showLegend is an example for generating legend in the figure of heatmap3 function. You can use your any plot functions to generate your own legend.
##' 
##' 
##' 
##' @inheritParams graphics::legend
##' @param lwd the line widths for lines appearing in the legend.
##' @param ... additional arguments passed on to \code{\link{legend}}
##' @export
##' @examples RowSideColors<-rep("steelblue2",nrow(mtcars))
##' RowSideColors[c(4:6,15:17,22:26,29)]<-"lightgoldenrod"
##' RowSideColors[c(1:3,19:21)]<-"brown1"
##' heatmap3(mtcars,scale="col",margins=c(2,10),RowSideColors=RowSideColors,legendfun=function() 
##' showLegend(legend=c("European","American","Japanese"),col=c("steelblue2","lightgoldenrod",
##' "brown1"),cex=1.5))
showLegend<-function(legend=c("Group A","Group B"),lwd=3,cex=1.1,col=c("red","blue"),...) {
	plot(0,xaxt="n",bty="n",yaxt="n",type="n",xlab="",ylab="")
	legend("topleft",legend=legend,lwd=lwd,col=col,bty="n",cex=cex,...)
}

##' showAnn
##' 
##' The function showAnn is an example for generating annotation figure in the result of heatmap3 function. You can use your any plot functions to generate your own annotation figure.
##' 
##' 
##' 
##' @param annData a data frame contains the annotation information for samples. It can only contain factor or numeric variables, and each row reprezent a sample with the same order of the columns in expression matrix.
##' @export
##' @examples annData<-data.frame(mtcars[,c("mpg","am","wt","gear")])
##' annData[,2]<-as.factor(annData[,2])
##' annData[,4]<-as.factor(annData[,4])
##' #Display annotation
##' \dontrun{
##' showAnn(annData)
##' }
##' #Heatmap with annotation
##' heatmap3(t(mtcars),ColSideAnn=annData,ColSideFun=function(x) 
##' showAnn(x),ColSideWidth=1.2,balanceColor=TRUE)
showAnn<-function(annData) {	
	#see how many lines were needed for factor annotation
	temp<-which(sapply(annData,class)=="factor")
	levelsAll<-sapply(annData[,temp,drop=F],levels)
	n<-length(unlist(levelsAll))
	totalLine<-(ncol(annData)-length(temp))*2+n
	plot(c(0-1/nrow(annData)/2,1+1/nrow(annData)/2),c(1,totalLine+1),type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",axes = FALSE,xaxs="i")
	
	offset<-1
	for (x in 1:ncol(annData)) {
		if (is.factor(annData[,x])) {
			temp<-factor2Matrix(annData[,x],colnames(annData)[x])
			image(seq(0,1,length.out=nrow(annData)),(1:ncol(temp))+offset-0.5,temp,col=c("white","black"),add=T)
			segments(seq(0-1/nrow(annData)/2,1+1/nrow(annData)/2,length.out=nrow(annData)+1),offset,seq(0-1/nrow(annData)/2,1+1/nrow(annData)/2,length.out=nrow(annData)+1),ncol(temp)+offset,col="white")
			abline(h=c(0,ncol(temp))+offset,col="white")
			mtext(2,at=(1:ncol(temp))+offset-0.5,text=colnames(temp),las=1)
			offset<-offset+ncol(temp)
		} else {
			points(seq(0,1,length.out=nrow(annData)),data2range(annData[,x],offset=offset))
			lines(seq(0,1,length.out=nrow(annData)),lowess(data2range(annData[,x]))$y+offset)
			mtext(2,at=offset+0.5,text=colnames(annData)[x],las=1)
			axis(4,at=range(data2range(annData[,x],offset=offset),na.rm=T),labels=round(range(annData[,x],na.rm=T),1),las=1)
			offset<-offset+2
		}
	}
	return(c(0.6,offset))
}

data2range<-function(x,minNew=0.1,maxNew=1.9,offset=0) {
	(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))*(maxNew-minNew)+minNew+offset
}

factor2Matrix<-function(factorData,colName) {
	factorMatrix<-matrix(0,nrow=length(factorData),ncol=length(levels(factorData)))
	temp<-cbind(1:nrow(factorMatrix),as.numeric(factorData))
	factorMatrix[temp]<-1
	colnames(factorMatrix)<-paste(colName,levels(factorData),sep="=")
	return(factorMatrix)
}
