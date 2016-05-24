PlotPathways<-function(Pathways,nRow=5,main=NULL,plottype="new",location=NULL){	
	
	if (!requireNamespace("MLP", quietly = TRUE)) {
		stop("MLP needed for this function to work. Please install it.",
				call. = FALSE)
	}
	
	if (!requireNamespace("biomaRt", quietly = TRUE)) {
		stop("biomaRt needed for this function to work. Please install it.",
				call. = FALSE)
	}
	
	if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
		stop("org.Hs.eg.db needed for this function to work. Please install it.",
				call. = FALSE)
	}
	
	
	plottypein<-function(plottype,location){
		if(plottype=="pdf" & !(is.null(location))){
			grDevices::pdf(paste(location,".pdf",sep=""))
		}
		if(plottype=="new"){
			grDevices::dev.new()
		}
		if(plottype=="sweave"){
			
		}
	}
	plottypeout<-function(plottype){
		if(plottype=="pdf"){
			grDevices::dev.off()
		}
	}
	#preparing data structure for the plotGOGraph
	colnames(Pathways)[3:ncol(Pathways)]=sub("mean_","",colnames(Pathways)[3:ncol(Pathways)])	
	#plot GOgraph
	plottypein(plottype,location)
	MLP::plotGOgraph(Pathways,nRow=nRow,main=main)
	plottypeout(plottype)
	
}


