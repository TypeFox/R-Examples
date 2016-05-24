tinGraphs <- function(res,DESIGN=NULL,x_axis=NULL,y_axis=NULL,inference.info=NULL,color.by.boots=TRUE,boot.cols=c('plum4','darkseagreen','firebrick3'),fi.col=NULL, fi.pch=NULL,fii.col=NULL,fii.pch=NULL,fj.col=NULL,fj.pch=NULL,col.offset=NULL,constraints=NULL,xlab=NULL,ylab=NULL,main=NULL,bootstrapBars=TRUE,correlationPlotter=TRUE,showHulls=0.95,biplots=FALSE){
		
	pca.types <- c('tepBADA')
	ca.types <- c('tepDICA')
	
	#A simple override/check. If someone puts in expoOutput class data, epGraphs will recognize it.
	if(class(res)[1] == "tinpoOutput"){
		if(length(res)==2){
			inference.info <- res$Inference.Data
		}
		res <- res$Fixed.Data
	}	
	
	#A simple override/check. If someone puts in texpoOutput class data, tepGraphs will recognize it.
	tepPlotInfo <- NULL
	if(class(res)[1] == "texpoOutput"){
		if(length(res)==2){
			tepPlotInfo <- res$Plotting.Data
		}
		res <- res$TExPosition.Data
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
	
	#perhaps make this stuff a function, or have TExPosition call all of tepGraphs.
	if(!(class(res)[1] %in% c(pca.types,ca.types))){
		stop("Unknown TExPosition class. Plotting has stopped.")
	}else{
		if(is.null(main)){
			main <- paste("Inferential ",deparse(substitute(res)),". Omni p=", inference.info$omni$p.val,sep="")
		}
		if(length(unlist(strsplit(main,"")))>40){
			main <- paste("Inferential Results. Omni p=", inference.info$omni$p.val,sep="")
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
		
		#tepPlotInfo check will look for proper colors & constraints, mostly.
		if(!is.null(tepPlotInfo)){
			if(!(nrow(res$fi)==nrow(tepPlotInfo$fi.col)) || !(nrow(res$fj)==nrow(tepPlotInfo$fj.col)) || !(nrow(res$fii)==nrow(tepPlotInfo$fii.col))){
				
				print('Dimension mismatch. tepPlotInfo will be reset, no hulls can be shown.')
				tepPlotInfo$fii.col <- NULL
				tepPlotInfo$fi.col <- NULL
				tepPlotInfo$fj.col <- NULL
				tepPlotInfo$constraints <- NULL
			}		
		}else{
			tepPlotInfo <- list(fii.col=NULL,fi.col=NULL,fj.col=NULL,constraints=NULL)
		}
		
		#fii.col, fi.col, fj.col, and constraints take precedence over tepPlotInfo. This is because epPlotInfo only exists via expoOutput.			
		if(is.null(fii.col) || is.null(fi.col) || nrow(fi.col)!=nrow(res$fi) || nrow(fii.col)!=nrow(res$fi)){
			if(is.null(tepPlotInfo$fii.col) || is.null(tepPlotInfo$fi.col)){
				if(is.null(DESIGN)){
					stop("fii.col and DESIGN are NULL. You must provide one or the other.")
				}else{
					#this will catch failures and stop.
					DESIGN <- texpoDesignCheck(DATA=NULL,DESIGN=DESIGN,make_design_nominal=FALSE)
					obs.cols <- createColorVectorsByDesign(DESIGN,offset=col.offset)
					fii.col <- obs.cols$oc
					fi.col <- obs.cols$gc
				}
			}else{
				fii.col <- tepPlotInfo$fii.col
				fi.col <- tepPlotInfo$fi.col
			}
		}
		if(is.null(fj.col) || nrow(fj.col)!=nrow(res$fj)){
			if(is.null(tepPlotInfo$fj.col)){
				fj.col <- createColorVectorsByDesign(matrix(1,nrow(res$fj),1),hsv=FALSE,offset=col.offset)$oc
			}else{
				fj.col <- tepPlotInfo$fj.col	
			}
		}
		
		if(is.null(fi.pch) || nrow(fi.pch)!=nrow(res$fi)){
			if(is.null(tepPlotInfo$fi.pch)){
				fi.pch <- as.matrix(rep(21,nrow(res$fi)))
			}else{
				fi.pch <- tepPlotInfo$fi.pch
			}
		}
		
		if(is.null(fii.pch) || nrow(fii.pch)!=nrow(res$fii)){
			if(is.null(tepPlotInfo$fii.pch)){
				fii.pch <- as.matrix(rep(21,nrow(res$fii)))
			}else{
				fii.pch <- tepPlotInfo$fii.pch
			}
		}
				
		if(is.null(fj.pch) || nrow(fj.pch)!=nrow(res$fj)){
			if(is.null(tepPlotInfo$fj.pch)){
				fj.pch <- as.matrix(rep(21,nrow(res$fj)))
			}else{
				fj.pch <- tepPlotInfo$fj.pch
			}
		}		
		
		if(is.null(constraints)){
			if(!is.null(tepPlotInfo$constraints)){
				constraints <- tepPlotInfo$constraints
			}
			#this is needed because if we switch axes, it could be different constraints.
			constraints <- calculateConstraints(results=res,x_axis=x_axis,y_axis=y_axis,constraints=constraints)			
		}		
		#by the time I get here, I should be guaranteed to have a fii.col, fi.col, fj.col, and constraints.
		
		fj.boot.tests <- rowSums(inference.info$boot.data$fj.boot.data$tests$sig.boot.ratios[,c(x_axis,y_axis)])
		fj.no.boot.axes <- which(fj.boot.tests == 0)
		fj.both.boot.axes <- which(fj.boot.tests == 2)
		fj.x.boot.axis <- which((inference.info$boot.data$fj.boot.data$tests$sig.boot.ratios[,c(x_axis)] - inference.info$boot.data$fj.boot.data$tests$sig.boot.ratios[,c(y_axis)])==1)
		fj.y.boot.axis <- which((inference.info$boot.data$fj.boot.data$tests$sig.boot.ratios[,c(y_axis)] - inference.info$boot.data$fj.boot.data$tests$sig.boot.ratios[,c(x_axis)])==1)
		
		
		fi.boot.tests <- rowSums(inference.info$boot.data$fi.boot.data$tests$sig.boot.ratios[,c(x_axis,y_axis)])
		fi.no.boot.axes <- which(fi.boot.tests == 0)
		fi.both.boot.axes <- which(fi.boot.tests == 2)
		fi.x.boot.axis <- which((inference.info$boot.data$fi.boot.data$tests$sig.boot.ratios[,c(x_axis)] - inference.info$boot.data$fi.boot.data$tests$sig.boot.ratios[,c(y_axis)])==1)
		fi.y.boot.axis <- which((inference.info$boot.data$fi.boot.data$tests$sig.boot.ratios[,c(y_axis)] - inference.info$boot.data$fi.boot.data$tests$sig.boot.ratios[,c(x_axis)])==1)
				
		orig.fi.col <- fi.col		
		#always do this -- user has no choice here.
		if(length(fj.no.boot.axes) != 0){
			fj.col[fj.no.boot.axes,1]  <- 'gray'
		}		
		if(length(fi.no.boot.axes) != 0){
			fi.col[fi.no.boot.axes,1]  <- 'gray'
		}				
		#everything here. -- DO I REALLY NEED ALL THESE?
		if(color.by.boots){
			if(is.null(boot.cols) || length(boot.cols) != 3){
				boot.cols <- c('plum4','darkseagreen','firebrick3')
			}				
			if(length(fj.x.boot.axis)!=0){
				fj.col[fj.x.boot.axis,1]  <- boot.cols[1]
			}
			if(length(fi.x.boot.axis)!=0){
				fi.col[fi.x.boot.axis,1]  <- boot.cols[1]
			}	
					
			if(length(fj.y.boot.axis)!=0){
				fj.col[fj.y.boot.axis,1]  <- boot.cols[2]
			}
			if(length(fi.y.boot.axis)!=0){
				fi.col[fi.y.boot.axis,1]  <- boot.cols[2]
			}			
			
			if(length(fj.both.boot.axes)!=0){
				fj.col[fj.both.boot.axes,1]  <- boot.cols[3]
			}
			if(length(fi.both.boot.axes)!=0){
				fi.col[fi.both.boot.axes,1]  <- boot.cols[3]
			}			
			
			fj.col.y <- fj.col.x <- fj.col	
			fj.col.y[fj.x.boot.axis,1] <- 'gray'
			fj.col.x[fj.y.boot.axis,1] <- 'gray'

			fi.col.y <- fi.col.x <- fi.col	
			fi.col.y[fi.x.boot.axis,1] <- 'gray'
			fi.col.x[fi.y.boot.axis,1] <- 'gray'			
		}		
		
		
		
		#OK SO WHAT IS NEXT?
		##a bit of fixed...
		fii.plot.info <- prettyPlot(res$fii,x_axis=x_axis,y_axis=y_axis,col=fii.col,pch=fii.pch,axes=TRUE,xlab=xlab,ylab=ylab,main=main,constraints=constraints,contributionCircles=FALSE,dev.new=TRUE)
		fi.plot.info <- prettyPlot(res$fi,x_axis=x_axis,y_axis=y_axis,col= orig.fi.col,pch=fi.pch,axes=FALSE,contributionCircles=TRUE,contributions=abs(inference.info$boot.data$fi.boot.data$tests$boot.ratios),dev.new=FALSE,new.plot=FALSE)
		if(showHulls > 0 && showHulls <= 1){
			colorDesign <- makeNominalData(fii.col)
			for(i in 1:nrow(res$fi)){
				peeledHull(res$fii[which(fii.col[, 1] == orig.fi.col[i,1]), ], x_axis = x_axis, y_axis = y_axis, percentage = showHulls, col = "black", lwd = 3)
				peeledHull(res$fii[which(fii.col[, 1] == orig.fi.col[i,1]), ], x_axis = x_axis, y_axis = y_axis, percentage = showHulls, col = orig.fi.col[i, ], lwd = 1)					
			}
		}
		
		
		#fii.plot.info <- prettyPlot(res$fii,x_axis=x_axis,y_axis=y_axis,col=fii.col,axes=TRUE,xlab=xlab,ylab=ylab,main=main,constraints=constraints,contributionCircles=FALSE,dev.new=TRUE)
		fi.plot.info <- prettyPlot(res$fi,x_axis=x_axis,y_axis=y_axis,col=fi.col,pch=fi.pch,axes=TRUE,xlab=xlab,ylab=ylab,main=main,constraints=constraints,contributionCircles=TRUE,contributions=abs(inference.info$boot.data$fi.boot.data$tests$boot.ratios),dev.new=TRUE)
		if(showHulls > 0 && showHulls <= 1){
			colorDesign <- makeNominalData(fii.col)
			for(i in 1:nrow(res$fi)){
				boot.items <- t(inference.info$boot.data$fi.boot.data$boots[i,,])
				peeledHull(boot.items, x_axis = x_axis, y_axis = y_axis, percentage = showHulls, col = "black", lwd = 5)
				peeledHull(boot.items, x_axis = x_axis, y_axis = y_axis, percentage = showHulls, col = fi.col[i,], lwd = 3)
			}
		}
		if(biplots){
			fj.plot.info <- prettyPlot(res$fj,x_axis=x_axis,y_axis=y_axis,col=fj.col,pch=fj.pch,axes=FALSE,contributionCircles=TRUE,contributions=abs(inference.info$boot.data$fj.boot.data$tests$boot.ratios),dev.new=FALSE,new.plot=FALSE)
		}else{
			fj.plot.info <- prettyPlot(res$fj,x_axis=x_axis,y_axis=y_axis,col=fj.col,pch=fj.pch,axes=TRUE,xlab=xlab,ylab=ylab,main=main,constraints=constraints,contributionCircles=TRUE,contributions=abs(inference.info$boot.data$fj.boot.data$tests$boot.ratios),dev.new=TRUE)	
		}
		
			
		if(bootstrapBars){
			# contributionBars(res$fi,res$ci,x_axis=x_axis,y_axis=y_axis,main=main,col=fi.plot.info$col)
			# contributionBars(res$fj,res$cj,x_axis=x_axis,y_axis=y_axis,main=main,col=fj.plot.info$col)
			prettyBars(inference.info$boot.data$fi.boot.data$tests$boot.ratios,axis=x_axis,fg.col=fi.col.x,dev.new=TRUE,threshold.line=TRUE,main=paste("Bootstrap Ratios Component: ",x_axis,sep=""),bg.lims=c(-inference.info$boot.data$fi.boot.data$tests$critical.value,inference.info$boot.data$fi.boot.data$tests$critical.value))
				
			prettyBars(inference.info$boot.data$fi.boot.data$tests$boot.ratios,axis=y_axis,fg.col=fi.col.y,dev.new=TRUE,horiz=FALSE,threshold.line=TRUE,main=paste("Bootstrap Ratios Component: ",y_axis,sep=""),bg.lims=c(-inference.info$boot.data$fi.boot.data$tests$critical.value,inference.info$boot.data$fi.boot.data$tests$critical.value))
			
			
			prettyBars(inference.info$boot.data$fj.boot.data$tests$boot.ratios,axis=x_axis,fg.col=fj.col.x,dev.new=TRUE,threshold.line=TRUE,main=paste("Bootstrap Ratios Component: ",x_axis,sep=""),bg.lims=c(-inference.info$boot.data$fj.boot.data$tests$critical.value,inference.info$boot.data$fj.boot.data$tests$critical.value))
				
			prettyBars(inference.info$boot.data$fj.boot.data$tests$boot.ratios,axis=y_axis,fg.col=fj.col.y,dev.new=TRUE,horiz=FALSE,threshold.line=TRUE,main=paste("Bootstrap Ratios Component: ",y_axis,sep=""),bg.lims=c(-inference.info$boot.data$fj.boot.data$tests$critical.value,inference.info$boot.data$fj.boot.data$tests$critical.value))					
			
		}		
		if(correlationPlotter && class(res)[1]%in%pca.types){
			correlationPlotter(res$X,res$fi,col=fj.col,x_axis=1,y_axis=2,xlab=xlab,ylab=ylab,main=main) 
		}											
		
	}	
}