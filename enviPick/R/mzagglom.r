mzagglom<-function(   MSlist,
                      dmzgap=10,
                      ppm=TRUE,
                      drtgap=500,
                      minpeak=4,
                      maxint=1E7,
                      progbar=FALSE
                      ){

	##############################################################################
	if(!length(MSlist)==8){stop("This is not an MSlist object")}
	if(!MSlist[[1]][[1]]){stop("MSlist empty or invalid. Use readMSdata to upload raw .mzML data first.")}
	if(!is.loaded("agglom")){stop(".Call to agglom failed; aborted.")}
	if(!is.loaded("indexed")){stop(".Call to indexed failed; aborted.")}
	if(!is.numeric(dmzgap)){stop("dmass must be numeric; aborted.")}
	if(!is.numeric(drtgap)){stop("dret must be numeric; aborted.")}
	if(!is.logical(ppm)){stop("ppm must be logical; aborted.")}
	############################################################################
	# reset if rerun ###########################################################
	MSlist[[5]]<-0;
	MSlist[[6]]<-0;
	MSlist[[7]]<-0;
	MSlist[[8]]<-0;
	MSlist[[4]][[2]][,5]<-rep(0,length(MSlist[[4]][[2]][,4]));
	MSlist[[4]][[2]][,6]<-rep(0,length(MSlist[[4]][[2]][,4]));
	MSlist[[4]][[2]][,7]<-rep(0,length(MSlist[[4]][[2]][,4]));
	##############################################################################
	# partition! #################################################################
	if(ppm){ppm2<-1}else{ppm2<-0};
	MSlist[[4]][[2]]<-MSlist[[4]][[2]][order(MSlist[[4]][[2]][,1],decreasing=FALSE),]
	if(progbar==TRUE){prog<-winProgressBar("Agglomerate...",min=0,max=3);
                      setWinProgressBar(prog, 0, title = "Agglomerate...", label = NULL);}
	part <- .Call("agglom",
			as.numeric(MSlist[[4]][[2]][,1]),
			as.numeric(MSlist[[4]][[2]][,3]),
			as.integer(ppm2),
			as.numeric(dmzgap),
			as.numeric(drtgap),
			PACKAGE="enviPick" 
		  )
	if(progbar==TRUE){setWinProgressBar(prog, 1, title = "Agglomerate...", label = NULL)}
	MSlist[[4]][[2]]<-MSlist[[4]][[2]][order(part,decreasing=FALSE),]
	part<-part[order(part,decreasing=FALSE)]
	maxit<-max(part)
	index <- .Call("indexed",
		as.integer(part),
		as.numeric(MSlist[[4]][[2]][,2]),
		as.integer(minpeak),
		as.numeric(maxint),
		as.integer(maxit),
		PACKAGE="enviPick" 
	)
	index<-index[index[,2]!=0,,drop=FALSE]
	colnames(index)<-c("start_ID","end_ID","number_peaks")
	MSlist[[5]]<-index
	partID<-.Call("partID",
		as.integer(index),
		as.integer(length(MSlist[[4]][[2]][,5])),
		PACKAGE="enviPick"  
	)
	MSlist[[4]][[2]][,5]<-partID
	if(progbar==TRUE){setWinProgressBar(prog, 2, title = "Agglomerate...", label = NULL)}
	##############################################################################
	MSlist[[3]][[1]]<-length(MSlist[[4]][[2]][,1]);
	MSlist[[3]][[2]]<-length(index[,3]);
	MSlist[[3]][[3]]<-sum(index[,3]);
	MSlist[[2]][[2]][1]<-as.character(dmzgap)
	MSlist[[2]][[2]][2]<-as.character(ppm)
	MSlist[[2]][[2]][3]<-as.character(drtgap)
	MSlist[[2]][[2]][4]<-as.character(minpeak)
	MSlist[[2]][[2]][5]<-as.character(maxint)  
	MSlist[[1]][[2]]<-TRUE;
	MSlist[[1]][[3]]<-FALSE;
	MSlist[[1]][[4]]<-FALSE;
	MSlist[[1]][[5]]<-FALSE;
	if(progbar==TRUE){close(prog);}
	##############################################################################
	return(MSlist)

}
