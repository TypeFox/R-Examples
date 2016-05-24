MultiTrainWithInventory<-function(model, CapacityHrs, ReserveHrs, RefillTime, DischargeCap=1, TurndownLimit=0.6, TurndownTime=1, ProgRpt=FALSE) {

	ShowProgress=FALSE					
		  OutputDF1=NULL				
		  OutputDF2=NULL				
		  ## since the first column of the input file is the Page entry we find the last one				
		  Pages=model[length(model[,1]),1]				
		  							
		    nextPageStarts<-1
## even proc.time has a problem on CRAN example run
if(ProgRpt==TRUE)  {			
	startTime<-proc.time()	
}			
		    				
	    for(p in 1:Pages)  {					
	    					
	    					
	    	thisPageStarts<-nextPageStarts				
	    	nextPageStarts=NULL				
	    	if(p<Pages) {				
	    		nextPage<-match(p+1,model[,1])			
	    		##nextPageStarts<-c(nextPageStarts,nextPage)			
	    		## not sure this variable was needed in this simpler case than OpDetail Dev required			
	    		nextPageStarts<-nextPage			
	    	}else{				
	    		## this makes the last page entry work			
	    		nextPage=length(model[,1])+1			
	    	}				
		    myMat<-as.matrix(model[thisPageStarts:(nextPage-1),4:ncol(model)])				
		    GenRate<-rowMeans(myMat)				
		    DurationVec<-model$Duration[thisPageStarts:(nextPage-1)]				
		    TimeVec<-model$Time[thisPageStarts:(nextPage-1)]				
		    NumTrains<-ncol(model)-3				
		Constraints<-c(CapacityHrs,ReserveHrs,RefillTime, DischargeCap,TurndownLimit,TurndownTime,NumTrains)				
		    				
		    ## this is the call to the C++ dll in the stosim library				
		      RcppList<-.Call("MultiTrainWithInventory",TimeVec, DurationVec, GenRate, Constraints, PACKAGE="stosim")				
		   				
		    lastline<-match(max(RcppList[[1]][,1]),RcppList[[1]][,1])				
		    RcppList[[1]]<-RcppList[[1]][1:lastline,]				
		    Page<-rep(p,length(RcppList[[1]][,1]))				
		    PageCol<-data.frame(Page)				
		    RcppList[[1]]<-cbind(PageCol,RcppList[[1]])				
		    OutputDF1<-rbind(OutputDF1,RcppList[[1]])				
		    				
		    Page<-rep(p,length(RcppList[[2]][,1]))				
		    PageCol<-data.frame(Page)				
		    RcppList[[2]]<-cbind(PageCol,RcppList[[2]])				
		    OutputDF2<-rbind(OutputDF2,RcppList[[2]])				
		    rm(RcppList)				
		    				
		    				
		    				
		    				
		    				
		    if(p==1)  { 
				if(ProgRpt==TRUE)  {			
				  oneCycle<-proc.time()				
				  TimeTest<-(oneCycle[3]-startTime[3])*Pages	  
				  ## still need to find the correct TimeTest value				
					  if(TimeTest > .5)  {				
					  pb <- tkProgressBar(title = "MultiTrainSingleBU Progress", min = 0,				
							   max = Pages, width = 300)				
					  ShowProgress=TRUE
					  }				
				}	
			}			  
		      if(ShowProgress==TRUE)  {				
		      setTkProgressBar(pb, p, label=paste( round(p/Pages*100, 0),"% done"))				
		    }				
		    				
		  ## return for next page				
		  }				
		  				
		  if(ShowProgress==TRUE)  {				
		  close(pb)				
		  }				
		  				
		  OutputList<-list(OutputDF1,OutputDF2)				
   return(OutputList)  						
}						
