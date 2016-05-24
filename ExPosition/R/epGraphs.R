epGraphs <-
function(res,x_axis=1,y_axis=2,epPlotInfo=NULL,DESIGN=NULL,fi.col=NULL,fi.pch=NULL,fj.col=NULL,fj.pch=NULL,col.offset=NULL,constraints=NULL,xlab=NULL,ylab=NULL,main=NULL,contributionPlots=TRUE,correlationPlotter=TRUE,biplots=FALSE,graphs=TRUE){

	pca.types <- c('epPCA','epMDS','epGPCA')
	ca.types <- c('epCA','epMCA')	
	
	#A simple override/check. If someone puts in expoOutput class data, epGraphs will recognize it.
	if(class(res)[1] == "expoOutput"){
		if(length(res)==2){
			epPlotInfo <- res$Plotting.Data
		}
		res <- res$ExPosition.Data
	}
	
	if(!(class(res)[1] %in% c(pca.types,ca.types))){
		stop("Unknown ExPosition class. Plotting has stopped.")
	}
	if(!is.null(epPlotInfo) && (class(epPlotInfo)[1] != "epGraphs")){
		stop("Unknown epPlotInfo class. Plotting has stopped.")
	}

	#epPlotInfo check will look for proper colors & constraints, mostly.
		###these NULL out because I'm going to test for them after this.
	if(!is.null(epPlotInfo)){
		if( !(nrow(res$fi)==nrow(epPlotInfo$fi.col)) ){
			print('$fi Dimension mismatch. epPlotInfo will be reset.')
			epPlotInfo <- list(fi.col=NULL,fi.pch=NULL,fj.col=NULL,fj.pch=NULL,constraints=NULL)
		}
		if( (!(class(res)[1]=='epMDS')) && !(nrow(res$fj)==nrow(epPlotInfo$fj.col)) ){
			print('$fj Dimension mismatch. epPlotInfo will be reset.')
			epPlotInfo <- list(fi.col=NULL,fi.pch=NULL,fj.col=NULL,fj.pch=NULL,constraints=NULL)
		}
	}else{
		epPlotInfo <- list(fi.col=NULL,fi.pch=NULL,fj.col=NULL,fj.pch=NULL,constraints=NULL)	
	}

	###Use this block to establish defaults.
	if(is.null(main)){
		main <- deparse(substitute(res))
	}		
	if(length(unlist(strsplit(main,"")))>40){
		main <- "Results"
	}		
	if(is.null(xlab)){
		xlab <- paste("Component ",x_axis," variance: ", round(res$t[x_axis],3), "%",sep="")
	}
	if(is.null(ylab)){
		ylab <- paste("Component ",y_axis," variance: ", round(res$t[y_axis],3), "%",sep="")
	}		
	if( (!is.null(col.offset)) && is.numeric(col.offset)){
		if(col.offset > 1){
			col.offset <- col.offset / as.numeric(paste(c(1,rep(0,nchar(as.character(col.offset)))),collapse=""))
		}
	}

	###establish fi.col
	if(length(fi.col)==1){
		fi.col <- as.matrix(rep(fi.col,nrow(res$fi)))
	}else if(is.null(fi.col) && !is.null(DESIGN)){
		fi.col <- createColorVectorsByDesign(DESIGN,offset=col.offset)$oc
	}else if(is.null(fi.col) && !is.null(epPlotInfo$fi.col)){
		fi.col <- epPlotInfo$fi.col 
	}
	###This is a final check.
	if(is.null(fi.col)){
		fi.col <- createColorVectorsByDesign(matrix(1,nrow(res$fi),1),offset=col.offset)$oc
	}
	if(nrow(fi.col)!=nrow(res$fi)){
		print('Incorrect fi.col. Creating default colors.')
		fi.col <- createColorVectorsByDesign(matrix(1,nrow(res$fi),1),offset=col.offset)$oc
	}
	
	###establish fi.pch	
	if(length(fi.pch)==1){
		fi.pch <- as.matrix(rep(fi.pch,nrow(res$fi)))
	}else if(is.null(fi.pch) && !is.null(epPlotInfo$fi.pch)){
		fi.pch <- epPlotInfo$fi.pch 
	}
	###This is a final check.
	if(is.null(fi.pch)){
		fi.pch <- as.matrix(rep(21,nrow(res$fi)))
	}
	if(nrow(fi.pch)!=nrow(res$fi)){
		print('Incorrect fi.pch. Creating default pch.')
		fi.pch <- as.matrix(rep(21,nrow(res$fi)))
	}
	
	#fjs
	if(class(res)[1]!='epMDS'){
		###establish fj.col
		if(length(fj.col)==1){
			fj.col <- as.matrix(rep(fj.col,nrow(res$fj)))
		}else if(is.null(fj.col) && !is.null(epPlotInfo$fj.col)){
			fj.col <- epPlotInfo$fj.col 
		}
		###This is a final check.
		if(is.null(fj.col)){
			fj.col <- createColorVectorsByDesign(matrix(1,nrow(res$fj),1),hsv=FALSE)$oc
		}
		if(nrow(fj.col)!=nrow(res$fj)){
			print('Incorrect fj.col. Creating default colors.')
			fj.col <- createColorVectorsByDesign(matrix(1,nrow(res$fj),1),hsv=FALSE)$oc
		}	
	
		###establish fj.pch	
		if(length(fj.pch)==1){
			fj.pch <- as.matrix(rep(fj.pch,nrow(res$fj)))
		}else if(is.null(fj.pch) && !is.null(epPlotInfo$fj.pch)){
			fj.pch <- epPlotInfo$fj.pch 
		}
		###This is a final check.
		if(is.null(fj.pch)){
			fj.pch <- as.matrix(rep(21,nrow(res$fj)))
		}
		if(nrow(fj.pch)!=nrow(res$fj)){
			print('Incorrect fj.pch. Creating default pch.')
			fj.pch <- as.matrix(rep(21,nrow(res$fj)))
		}
	}	
	
	###establish constraints
	if(is.null(constraints) && !is.null(epPlotInfo$constraints)){
		constraints <- epPlotInfo$constraints
	}
	##constraints always need to be checked.
	constraints <- calculateConstraints(results=res,x_axis=x_axis,y_axis=y_axis,constraints=constraints)
	
	if(graphs){
		fi.plot.info <- prettyPlot(res$fi,x_axis=x_axis,y_axis=y_axis,col=fi.col,axes=TRUE,xlab=xlab,ylab=ylab,main=main,constraints=constraints,pch=fi.pch,contributionCircles=TRUE,contributions=res$ci,dev.new=TRUE)
		if(!(class(res)[1]=='epMDS')){
			if(biplots){
				fj.plot.info <- prettyPlot(res$fj,x_axis=x_axis,y_axis=y_axis,col=fj.col,axes=FALSE,contributionCircles=TRUE,contributions=res$cj,pch=fj.pch,dev.new=FALSE,new.plot=FALSE)
			}else{
				fj.plot.info <- prettyPlot(res$fj,x_axis=x_axis,y_axis=y_axis,col=fj.col,axes=TRUE,xlab=xlab,ylab=ylab,main=main,constraints=constraints,pch=fj.pch,contributionCircles=TRUE,contributions=res$cj,dev.new=TRUE)		
			}	
		}
		if(contributionPlots){
			contributionBars(res$fi,res$ci,x_axis=x_axis,y_axis=y_axis,main=main,col=fi.plot.info$col)
			if(!(class(res)[1]=='epMDS')){
				contributionBars(res$fj,res$cj,x_axis=x_axis,y_axis=y_axis,main=main,col=fj.plot.info$col)
			}
		}		
		if(correlationPlotter && class(res)[1]%in%pca.types){
			if(class(res)[1]=='epMDS'){
				correlationPlotter(res$X,res$fi,col=fi.col,pch=fi.pch,x_axis=x_axis,y_axis=y_axis,xlab=xlab,ylab=ylab,main=main) 
			}else{
				correlationPlotter(res$X,res$fi,col=fj.col,pch=fj.pch,x_axis=x_axis,y_axis=y_axis,xlab=xlab,ylab=ylab,main=main) 
			}
		}
	}								
	
	#this happens whether I graph, or not. 
	if(class(res)[1]=='epMDS'){
		epPlotInfo <- list(fi.col=fi.col,fi.pch=fi.pch,fj.col=fi.col,fj.pch=fi.pch,constraints=constraints)	
	}else{
		epPlotInfo <- list(fi.col=fi.col,fi.pch=fi.pch,fj.col=fj.col,fj.pch=fj.pch,constraints=constraints)
	}

	class(epPlotInfo) <- c("epGraphs", "list")
	return(epPlotInfo)	
}