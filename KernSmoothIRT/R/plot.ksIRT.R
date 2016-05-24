plot.ksIRT <-
function(x, plottype = c("OCC","EIS","density","expected","sd","triangle","tetrahedron","RCC","EISDIF","OCCDIF","PCA","expectedDIF","densityDIF"), items= "all", subjects, axistype = c("scores","distribution"), alpha, main, xlab, ylab, xlim, ylim, cex,...){

	oldops <- options()
	oldpar <- par(no.readonly = TRUE)
	

	plottype <- match.arg(arg = plottype, choices = c("OCC","EIS","density","expected","sd","triangle","tetrahedron","RCC","EISDIF","OCCDIF","PCA","expectedDIF","densityDIF"))
	axistype <- match.arg(arg=axistype, choices = c("scores","distribution"))
	
	if(items[1]=="all"){items<-1:x$nitem}
	if(missing(xlab)){xlab0<-"Scores"}
	if(plottype=="expected"){axistype="distribution"}
	if(plottype=="expectedDIF"){axistype="scores"}

	if(axistype=='distribution'){

		axis<-x$evalpoints
		quants<-x$subjthetasummary
		if(missing(xlab)){
			xlab<-paste(simpleCap(x$thetadist[[1]]), paste(unlist(x$thetadist[-1]),collapse=" "), "Quantiles",sep=" ")
		}
	}
	else{

		axis<-x$expectedscores
		quants<-x$subjscoresummary
		if(missing(xlab)){	
			xlab<-"Expected Score"
		}
	}

	switch(plottype,
			density = densityplot(x,xlim,ylim,xlab,ylab,main,...),
			EIS = ICCplot(x,items,alpha,axis=axis,quants=quants,main,xlab,ylab,xlim,ylim,cex,...),
			OCC = OCCplot(x,items,alpha,axis=axis,quants=quants,main,xlab,ylab,xlim,ylim,...),
			expected = expectedplot(x,axis=axis,quants=x$subjthetasummary, main, xlim,ylim, xlab, ylab, ...),
			sd = sdplot(x,axis=axis,quants=quants,main, xlab, ylab=ylab,...),
			triangle = Simplexplot(x,items,main, ...),
			tetrahedron = Tetplot(x,items,main, ...),
			RCC = Credplot(x,axis=axis,quants=quants,xlab=xlab,subjects=subjects,xlim,...),
			EISDIF = ICCDIFplot(x,items,alpha,axis=axis,quants=quants,main,xlab,ylab,xlim,ylim,cex,...),
			OCCDIF = OCCDIFplot(x,items,alpha,axis=axis,quants=quants,main,xlab,ylab,xlim,ylim,...),
			PCA = PCAplot(x, ...),
			expectedDIF = expectedDIFplot(x,axis=axis,quants=quants, main, xlim,ylim, xlab, ylab, ...),
			densityDIF = densityDIFplot(x,xlim,ylim,xlab=xlab0,ylab,main,...),
			print("plottype not recognized")
	)
	
	options(oldops)
	par(oldpar)
	
	
}

