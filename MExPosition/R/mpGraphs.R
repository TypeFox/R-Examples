mpGraphs <- function(res, table, DESIGN = NULL, x_axis=1, y_axis=2, fi.col=NULL, fj.col=NULL, table.col = NULL, col.offset = NULL, constraints=NULL, xlab = NULL, ylab = NULL, main = NULL, graphs=TRUE) 
{	mpPlotInfo <- NULL
	if(class(res)[1] == "mexpoOutput")
	{	if(length(res) == 2)
		{	mpPlotInfo <- res$Plotting.Data
		}
		res <- res$mexPosition.Data
	} 

	if(is.null(main))
	{	main <- deparse(substitute(res))
		if(length(unlist(strsplit(main,"")))>40)
		{	main <- "Results"	
		}
	}
	
	if(!(class(res)[1] %in% c('mpSTATIS','mpDISTATIS','mpMFA','mpCANOSTATIS','mpANISOSTATIS','mpCOVSTATIS','mpPTA','mpKPlus1STATIS','mpDOACT.STATIS')))
	{	stop('Unknown MExPosition class. Plotting has stopped.')
		
	} 

	if(!(class(res)[1] %in% c('mpDOACT.STATIS')))
	{	if(is.null(xlab))
		{	xlab.innerproduct <- paste("Component ",x_axis, " variance:", round(res$InnerProduct$t[x_axis],3),"%",sep="")
			xlab.table <- paste("Component ",x_axis, " variance:", round(res$Table$t[x_axis],3),"%",sep="")
		}
		if(is.null(ylab))
		{	ylab.innerproduct <- paste("Component ",y_axis, " variance:", round(res$InnerProduct$t[y_axis],3),"%",sep="")
		ylab.table <- paste("Component ",y_axis, " variance:", round(res$Table$t[y_axis],3),"%",sep="")
		}
	}
	if((class(res)[1] %in% c('mpDOACT.STATIS')))
	{	if(is.null(xlab))
		{	xlab.innerproduct <- paste("Component ",x_axis, " variance:", round(res$InnerProduct$t[x_axis],3),"%",sep="")
			xlab.table.1 <- paste("Component ",x_axis, " variance:", round(res$Table$t.1[x_axis],3),"%",sep="")
			xlab.table.2 <- paste("Component ",x_axis, " variance:", round(res$Table$t.2[x_axis],3),"%",sep="")
		}
		if(is.null(ylab))
		{	ylab.innerproduct <- paste("Component ",y_axis, " variance:", round(res$InnerProduct$t[y_axis],3),"%",sep="")
			ylab.table.1 <- paste("Component ",y_axis, " variance:", round(res$Table$t.1[y_axis],3),"%",sep="")
			ylab.table.2 <- paste("Component ",y_axis, " variance:", round(res$Table$t.2[y_axis],3),"%",sep="")
		}
	}



	if(!is.null(mpPlotInfo))
	{	if(!(nrow(res$fi)==nrow(mpPlotInfo$ficol)) || !(nrow(res$fj)==nrow(mpPlotInfo$fj.col)))
		{	print('Dimensions mismatch. mpPlotInfo will be reset')
			mpPlotInfo$fi.col <- NULL
			mpPlotInfo$fj.col <- NULL
			mpPlotInfo$table.col <- NULL
			mpPlotInfo$constraints <- NULL
		}
	}
	else
	{	mpPlotInfo <- list(fi.col = NULL, fj.col = NULL, table.col = NULL, constraints = NULL)
	}

	if(is.null(fi.col) || nrow(fi.col) != nrow(res$fi))
	{	if(is.null(mpPlotInfo$fi.col))
		{	if(is.null(DESIGN))
			{	fi.col <- createColorVectorsByDesign(matrix(1,nrow(res$InnerProduct$fi),1),offset=col.offset)$oc
				fj.col <- createColorVectorsByDesign(t(table),offset=col.offset)$oc
				table.col <- createColorVectorsByDesign(t(table),offset=col.offset)$gc
			}
			else
			{	fi.col <- createColorVectorsByDesign(DESIGN, offset=col.offset)$oc
				fj.col <- createColorVectorsByDesign(t(table), offset=col.offset)$oc
				table.col <- createColorVectorsByDesign(t(table), offset=col.offset)$gc
			}
		} 
	
		else
		{	fi.col <- mpPlotInfo$fi.col
			fj.col <- mpPlotInfo$fj.col
			table.col <- mpPlotInfo$table.col
		}
	
	}

	if(is.null(constraints))
	{	if(!is.null(mpPlotInfo$constraints))
		{	constraints <- mpPlotInfo$constraints
		}
		constraints <- calculateConstraints(results=res$InnerProduct, x_axis = x_axis, y_axis=y_axis, constraints = constraints)
	}
	
	if(graphs)
	{	if((class(res)[1] %in% c('mpSTATIS','mpMFA','mpCANOSTATIS','mpANISOSTATIS','mpPTA','mpKPlus1STATIS')))	
		{	innerproduct.fi.plot.info <- prettyPlot(res$InnerProduct$fi, x_axis=x_axis, y_axis=y_axis, col=table.col, xlab=xlab.innerproduct, ylab=ylab.innerproduct, 
												main = paste("Inner Product",main), contributionCircles=TRUE, contributions = res$InnerProduct$ci, dev.new=TRUE)
		
			compromise.fi.plot.info <- prettyPlot(res$Compromise$compromise.fi, x_axis=x_axis, y_axis=y_axis, col=fi.col, xlab=xlab.table, ylab=ylab.table, 
												main = paste("Compromise",main), contributionCircles=TRUE, contributions = res$Compromise$compromise.ci,dev.new=TRUE)
		
			partial.plot.info <- prettyPlot(res$Table$partial.fi, x_axis=x_axis, y_axis=y_axis, col=as.matrix(rep(fi.col,res$Overview$num.groups)), 
												xlab=xlab.table, ylab=ylab.table, main = paste("Partial Scores",main),dev.new=TRUE)	

			fi.plot.info <- prettyPlot(res$Compromise$compromise.fi, x_axis=x_axis, y_axis=y_axis, col=fi.col,contributionCircles=TRUE, contributions = res$Compromise$compromise.ci, dev.new=FALSE, new.plot=FALSE)												
																	
			loading.plot.info <- prettyPlot(res$Table$Q, x_axis=x_axis, y_axis=y_axis, col=fj.col, xlab=xlab.table, ylab=ylab.table, 
												main = paste("Loadings",main),contributions = res$Table$cj, contributionCircles=TRUE, dev.new=TRUE)
		}
		
		if(class(res)[1] %in% c('mpDISTATIS','mpCOVSTATIS'))
		{	innerproduct.fi.plot.info <- prettyPlot(res$InnerProduct$fi, x_axis=x_axis, y_axis=y_axis, col=table.col, xlab=xlab.innerproduct, ylab=ylab.innerproduct, 
												main = paste("Inner Product",main), contributionCircles=TRUE, contributions = res$InnerProduct$ci, dev.new=TRUE)
		
			compromise.fi.plot.info <- prettyPlot(res$Compromise$compromise.fi, x_axis=x_axis, y_axis=y_axis, col=fi.col, xlab=xlab.table, ylab=ylab.table, 
												main = paste("Compromise",main), contributionCircles=TRUE, contributions = res$Compromise$compromise.ci, dev.new=TRUE)
				
												
			partial.plot.info <- prettyPlot(res$Table$partial.fi, x_axis=x_axis, y_axis=y_axis, col=as.matrix(rep(fi.col,res$Overview$num.groups)), 
			 									xlab=xlab.table, ylab=ylab.table, main = paste("Partial Scores",main), dev.new=TRUE)

			fi.plot.info <- prettyPlot(res$Compromise$compromise.fi, x_axis=x_axis, y_axis=y_axis, col=fi.col, contributionCircles=TRUE, contributions = res$Compromise$compromise.ci, dev.new=FALSE, new.plot=FALSE)
		
			loading.plot.info <- prettyPlot(res$Table$Q, x_axis=x_axis, y_axis=y_axis, col=fj.col, xlab=xlab.table, ylab=ylab.table, main = paste("Loadings",main),contributions = res$Table$cj, contributionCircles=TRUE, dev.new=TRUE)
		}

		if(class(res)[1] %in% c('mpDOACT.STATIS'))
		{	innerproduct.fi.plot.info <- prettyPlot(res$InnerProduct$fi, x_axis=x_axis, y_axis=y_axis, col=table.col, xlab=xlab.innerproduct, ylab=ylab.innerproduct, 
												main = paste("Inner Product",main), contributionCircles=TRUE, contributions = res$InnerProduct$ci, dev.new=TRUE)
		
			compromise.fi.1.plot.info <- prettyPlot(res$Compromise$compromise.fi.1, x_axis=x_axis, y_axis=y_axis, col=fi.col, xlab=xlab.table.1, ylab=ylab.table.1, 
												main = paste("Compromise (Dataset 1)",main), contributionCircles=TRUE, contributions = res$Compromise$compromise.ci.1, dev.new=TRUE)

			compromise.fi.2.plot.info <- prettyPlot(res$Compromise$compromise.fi.2, x_axis=x_axis, y_axis=y_axis, col=fi.col, xlab=xlab.table.2, ylab=ylab.table.2, 
												main = paste("Compromise (Dataset 2)",main), contributionCircles=TRUE, contributions = res$Compromise$compromise.ci.2, dev.new=TRUE)
		
			partial.1.plot.info <- prettyPlot(res$Table$partial.fi.1, x_axis=x_axis, y_axis=y_axis, col=as.matrix(rep(fi.col,res$Overview$num.groups.1)), 
												xlab=xlab.table.1, ylab=ylab.table.1, main = paste("Partial Scores (Dataset 1)",main), dev.new=TRUE)	

			fi.1.plot.info <- prettyPlot(res$Compromise$compromise.fi.1, x_axis=x_axis, y_axis=y_axis, col=fi.col,contributionCircles=TRUE, contributions = res$Compromise$compromise.ci.1, dev.new=FALSE, new.plot=FALSE)												
			
			partial.2.plot.info <- prettyPlot(res$Table$partial.fi.2, x_axis=x_axis, y_axis=y_axis, col=as.matrix(rep(fi.col,res$Overview$num.groups.2)), 
												xlab=xlab.table.2, ylab=ylab.table.2, main = paste("Partial Scores (Dataset 2)",main), dev.new=TRUE)	

			fi.1.plot.info <- prettyPlot(res$Compromise$compromise.fi.2, x_axis=x_axis, y_axis=y_axis, col=fi.col,contributionCircles=TRUE, contributions = res$Compromise$compromise.ci.2, dev.new=FALSE, new.plot=FALSE)												
			
			loading.1.plot.info <- prettyPlot(res$Table$Q.1, x_axis=x_axis, y_axis=y_axis, col=fj.col, xlab=xlab.table.1, ylab=ylab.table.1, 
												main = paste("Loadings (Dataset 1)",main),contributions = res$Table$cj.1, contributionCircles=TRUE, dev.new=TRUE)

			loading.2.plot.info <- prettyPlot(res$Table$Q.2, x_axis=x_axis, y_axis=y_axis, col=fj.col, xlab=xlab.table.1, ylab=ylab.table.1, 
												main = paste("Loadings (Dataset 1)",main),contributions = res$Table$cj.2, contributionCircles=TRUE, dev.new=TRUE)

		}
	}

		
	mpPlotInfo <- list(fi.col = fi.col, fj.col = fj.col, table.col=table.col, constraints = constraints)
	class(mpPlotInfo) <- c("mpGraphs","list")
	return(mpPlotInfo)
}
