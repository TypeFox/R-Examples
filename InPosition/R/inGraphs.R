inGraphs <- function(res,DESIGN=NULL,x_axis=NULL,y_axis=NULL,inference.info=NULL,color.by.boots=TRUE,boot.cols=c('plum4','darkseagreen','firebrick3'), fi.col=NULL,fi.pch=NULL,fj.col=NULL,fj.pch=NULL,col.offset=NULL,constraints=NULL,xlab=NULL,ylab=NULL,main=NULL,bootstrapBars=TRUE,correlationPlotter=TRUE,biplots=FALSE){

##update this to get the pchs.

	pca.types <- c('epPCA','epMDS','epGPCA')
	ca.types <- c('epCA','epMCA')

	#A simple override/check. If someone puts in expoOutput class data, epGraphs will recognize it.
	if(class(res)[1] == "inpoOutput"){
		if(length(res)==2){
			inference.info <- res$Inference.Data
		}
		res <- res$Fixed.Data
	}

	#A simple override/check. If someone puts in expoOutput class data, epGraphs will recognize it.
	epPlotInfo <- NULL
	if(class(res)[1] == "expoOutput"){
		if(length(res)==2){
			epPlotInfo <- res$Plotting.Data
		}
		res <- res$ExPosition.Data
	}		
	
	if(is.null(inference.info)){
		stop("inGraphs requires inference.info")
	}##could use some more checks...	

	component.p.order <- order(inference.info$components$p.vals)
	which.axes <- sort(component.p.order[1:2])	
	
	if(is.null(x_axis)){
		x_axis <- which.axes[1]
	}
	if(is.null(y_axis)){
		y_axis <- which.axes[2]		
	}
		
	#if the override/check fails, it skips to this. Which means it was internal to ep*()
	if(!(class(res)[1] %in% c(pca.types,ca.types))){
		stop("Unknown ExPosition class. Plotting has stopped.")
	}else{
		###Use this block to establish defaults.
		if(is.null(main)){
			main <- deparse(substitute(res))
			main <- paste("Inferential",main,sep=" ")
		}		
		if(length(unlist(strsplit(main,"")))>40){
			main <- "Inferential Results"
		}
		if(class(res)[1] %in% ca.types){
			main <- paste(main," Omni p=", inference.info$omni$p.val,sep="")
		}	
		if(is.null(xlab)){
			xlab <- paste("Component ",x_axis," variance: ", round(res$t[x_axis],3), "%, p=", inference.info$components$p.vals[x_axis],sep="")
		}else{
			xlab <- paste(xlab,"; p=", inference.info$components$p.vals[x_axis],sep="")			
		}
		if(is.null(ylab)){
			ylab <- paste("Component ",y_axis," variance: ", round(res$t[y_axis],3), "%, p=", inference.info$components$p.vals[y_axis],sep="")
		}else{
			ylab <- paste(ylab,"; p=", inference.info$components$p.vals[y_axis],sep="")
		}
		#epPlotInfo check will look for proper colors & constraints, mostly.
		if(!is.null(epPlotInfo)){
			if(class(res)[1]=='epMDS'){
				if(!(nrow(res$fi)==nrow(epPlotInfo$fi.col))){
					print('Dimension mismatch. epPlotInfo will be reset.')
					epPlotInfo$fi.col <- NULL
					epPlotInfo$fj.col <- NULL
					epPlotInfo$constraints <- NULL
				}
			}else{
				if(!(nrow(res$fi)==nrow(epPlotInfo$fi.col)) || !(nrow(res$fj)==nrow(epPlotInfo$fj.col))){
					print('Dimension mismatch. epPlotInfo will be reset.')
					epPlotInfo$fi.col <- NULL
					epPlotInfo$fj.col <- NULL
					epPlotInfo$constraints <- NULL
				}
			}
		}else{
			epPlotInfo <- list(fi.col=NULL,fj.col=NULL,constraints=NULL)	
		}

		#fi.col, fj.col, and constraints take precedence over epPlotInfo. This is because epPlotInfo only exists via expoOutput.					
		if(is.null(fi.col) || nrow(fi.col)!=nrow(res$fi)){
			if(is.null(epPlotInfo$fi.col)){
				if(is.null(DESIGN)){
					fi.col <- createColorVectorsByDesign(matrix(1,nrow(res$fi),1),offset=col.offset)$oc
				}else{
					fi.col <- createColorVectorsByDesign(DESIGN,offset=col.offset)$oc
				}
			}else{
				fi.col <- epPlotInfo$fi.col
			}
		}

		if(is.null(fi.pch) || nrow(fi.pch)!=nrow(res$fi)){
			if(is.null(epPlotInfo$fi.pch)){
				fi.pch <- fi.pch <- as.matrix(rep(21,nrow(res$fi)))
			}else{
				fi.pch <- epPlotInfo$fi.pch
			}
		}		
		
		if(class(res)[1]!='epMDS'){
			if(is.null(fj.col) || nrow(fj.col)!=nrow(res$fj)){
				if(is.null(epPlotInfo$fj.col)){
					fj.col <- createColorVectorsByDesign(matrix(1,nrow(res$fj),1),hsv=FALSE,offset=col.offset)$oc
				}else{
					fj.col <- epPlotInfo$fj.col	
				}
			}
			if(is.null(fj.pch) || nrow(fj.pch)!=nrow(res$fj)){
				if(is.null(epPlotInfo$fi.pch)){
					fj.pch <- fj.pch <- as.matrix(rep(21,nrow(res$fj)))
				}else{
					fj.pch <- epPlotInfo$fj.pch
				}
			}			
		}		
		
		if(is.null(constraints)){
			if(!is.null(epPlotInfo$constraints)){
				constraints <- epPlotInfo$constraints
			}
			#this is needed because if we switch axes, it could be different constraints.
			constraints <- calculateConstraints(results=res,x_axis=x_axis,y_axis=y_axis,constraints=constraints)			
		}
		#by the time I get here, I should be guaranteed to have a fi.col, fj.col, and constraints.		
		#this function will always graph.
		
		
		#fj.boot.cols <- fj.col
		boot.tests <- rowSums(inference.info$fj.boots$tests$sig.boot.ratios[,c(x_axis,y_axis)])
		no.boot.axes <- which(boot.tests == 0)
		both.boot.axes <- which(boot.tests == 2)
		#x.boot.axis <- which(inference.info$fj.boots$sig.boot.ratios[,c(x_axis)]==1)
		x.boot.axis <- which((inference.info$fj.boots$tests$sig.boot.ratios[,c(x_axis)] - inference.info$fj.boots$tests$sig.boot.ratios[,c(y_axis)])==1)
		#y.boot.axis <- which(inference.info$fj.boots$sig.boot.ratios[,c(y_axis)]==1)
		y.boot.axis <- which((inference.info$fj.boots$tests$sig.boot.ratios[,c(y_axis)] - inference.info$fj.boots$tests$sig.boot.ratios[,c(x_axis)])==1)
	
		#always do this -- user has no choice here.
		if(length(no.boot.axes) != 0){
			fj.col[no.boot.axes,1]  <- 'gray'
		}		
		#everything here.
		if(color.by.boots){
			if(is.null(boot.cols) || length(boot.cols) != 3){
				boot.cols <- c('plum4','darkseagreen','firebrick3')
			}				
			if(length(x.boot.axis)!=0){
				fj.col[x.boot.axis,1]  <- boot.cols[1]
			}
			if(length(y.boot.axis)!=0){
				fj.col[y.boot.axis,1]  <- boot.cols[2]
			}
			if(length(both.boot.axes)!=0){
				fj.col[both.boot.axes,1]  <- boot.cols[3]
			}
			fj.col.y <- fj.col.x <- fj.col			
			fj.col.y[x.boot.axis,1] <- 'gray'
			fj.col.x[y.boot.axis,1] <- 'gray'			
		}
		
		fi.plot.info <- prettyPlot(res$fi,x_axis=x_axis,y_axis=y_axis,col=fi.col,axes=TRUE,xlab=xlab,ylab=ylab,main=main,constraints=constraints,pch=fi.pch,contributionCircles=TRUE,contributions=res$ci,dev.new=TRUE)
		
		if(!(class(res)[1]=='epMDS')){
			if(biplots){
				fj.plot.info <- prettyPlot(res$fj,x_axis=x_axis,y_axis=y_axis,col=fj.col,axes=FALSE,pch=fj.pch,contributionCircles=TRUE,contributions=abs(inference.info$fj.boots$tests$boot.ratios),dev.new=FALSE,new.plot=FALSE)
			}else{
				fj.plot.info <- prettyPlot(res$fj,x_axis=x_axis,y_axis=y_axis,col=fj.col,axes=TRUE,xlab=xlab,ylab=ylab,main=main,constraints=constraints,pch=fj.pch,contributionCircles=TRUE,contributions=abs(inference.info$fj.boots$tests$boot.ratios),dev.new=TRUE)		
			}
		}
		if(bootstrapBars){
			if(!(class(res)[1]=='epMDS')){
				prettyBars(inference.info$fj.boots$tests$boot.ratios,axis=x_axis,fg.col=fj.col.x,dev.new=TRUE,threshold.line=TRUE,main=paste("Bootstrap Ratios Component: ",x_axis,sep=""),bg.lims=c(-inference.info$fj.boots$tests$critical.value,inference.info$fj.boots$tests$critical.value))
				
				prettyBars(inference.info$fj.boots$tests$boot.ratios,axis=y_axis,fg.col=fj.col.y,dev.new=TRUE,horiz=FALSE,threshold.line=TRUE,main=paste("Bootstrap Ratios Component: ",y_axis,sep=""),bg.lims=c(-inference.info$fj.boots$tests$critical.value,inference.info$fj.boots$tests$critical.value))
				
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
	
}