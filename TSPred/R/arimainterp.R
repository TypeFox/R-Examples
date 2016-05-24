arimainterp <-
function(TimeSeries,n.ahead,extrap=TRUE,xreg=NULL,newxreg=NULL,se.fit=FALSE){
    if(is.null(TimeSeries) | ncol(TimeSeries)<2) stop("TimeSeries is required and must have positive length with at least 2 columns")
    if(is.null(n.ahead)) stop("The number of values to be predicted is unknown")
    
	FrstBlocksCols <- c(1:(ncol(TimeSeries)-1))
    LstBlocksCols <- c(2:(ncol(TimeSeries)))
    
    FrstBlocks <- data.frame(TimeSeries[,FrstBlocksCols])
	RevLstBlocks <- data.frame(TimeSeries[order(nrow(TimeSeries):1),LstBlocksCols])
	
	N1 <- ceiling(n.ahead/2)
	N2 <- floor(n.ahead/2)
    
	namescol <- NULL
	for (i in FrstBlocksCols) namescol[i] <- paste('Block',i)	
	
	reg <- xreg
	newreg <- newxreg
	if(!is.null(xreg)) reg <- data.frame(xreg[FrstBlocksCols])
	if(!is.null(newxreg)) newreg <- mapply(head,data.frame(newxreg[FrstBlocksCols]),N1)
	FrstBlocksPredictions <- marimapred(FrstBlocks,n.ahead=N1,na.action=na.omit,xreg=reg,newxreg=newreg,se.fit=se.fit)
	
	if(!is.null(xreg)) reg <- sapply(data.frame(xreg[LstBlocksCols]),rev)
	if(!is.null(newxreg)) newreg <- apply( data.frame(mapply(tail,data.frame(newxreg[FrstBlocksCols]),N2)) ,2,rev)
	RevFrstBlocksPredictions <- marimapred(RevLstBlocks,n.ahead=N2,na.action=na.omit,xreg=reg,newxreg=newreg,se.fit=se.fit)
	
	RevFrstBlocksPredictions <- apply(RevFrstBlocksPredictions,2,rev)

    InterpPredictions <- rbind(FrstBlocksPredictions,RevFrstBlocksPredictions)
    colnames(InterpPredictions) <- namescol
	
    #Extrapolation phase
	if(extrap){		
		LstBlock <- TimeSeries[ncol(TimeSeries)]
		
		reg <- xreg
		newreg <- newxreg
		if(!is.null(xreg)) reg <- xreg[ncol(TimeSeries)]
		if(!is.null(newxreg)) newreg <- newxreg[ncol(TimeSeries)]

		LstBlockExtrapPredictions <- marimapred(LstBlock,n.ahead=n.ahead,na.action=na.omit,xreg=reg,newxreg=newreg,se.fit=se.fit)

		InterpPredictions <- cbind(InterpPredictions,LstBlockExtrapPredictions)
		colnames(InterpPredictions) <- c(namescol, paste('Block',length(namescol)+1))
	}

    return (InterpPredictions)
}