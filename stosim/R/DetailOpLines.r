				
 ## DetailOplines.r file				
 ##				
 ## Author: Jacob Ormerod				
 ##        ( c ) 2011-2014 OpenReliability.org				
##				
				
 DetailOpLines<-function(Model,Names=NULL, ProgRpt=FALSE)  {				
				
	## Initialize a ShowProgress boolean object			
	ShowProgress=FALSE			
	## Initialize the output dataframe			
	OutputDF=NULL			
	if(length(Names)>0) {			
		## assure that Names vector is same length as Model			
		if(length(Names) != length(Model) )  {			
			stop("Error in stosim:  Names count does not equal number of OpLines in Model")		
		}			
	}else{
		for(x in 1:length(Model))  {
			Names<-c(Names,paste("Train",as.character(x)))
			}
	}
	
	## error to assure input SimHistories are same size in terms of pages and hours per page			
	TestDF<-Model[[1]]			
	Pages<-TestDF$Page[length(TestDF$Page)]			
	HoursPerPage<-TestDF$Time[length(TestDF$Time)]			
				
	if(length(Model)>1)   {			
	for(sh in 2:length(Model))  {			
		TestDF<-Model[[sh]]		
		if(TestDF$Page[length(TestDF$Page)] != Pages)  {		
			stop("Error in stosim:  Simulation Histories are not consistent")	
		}		
		if(TestDF$Time[length(TestDF$Time)] != HoursPerPage)  {		
			stop("Error in stosim:  Simulation Histories are not consistent")	
		}		
	}			
	}			
				
	nextPageStarts<-rep(1,length(Model))			
				
## even proc.time has a problem on CRAN example run
if(ProgRpt==TRUE)  {			
	startTime<-proc.time()	
}				
				
for(p in 1:Pages)  {				
	thisPageStarts<-nextPageStarts			
	nextPageStarts=NULL			
				
				
	LengthsVec=NULL			
	DurationsVec=NULL			
	TimesVec=NULL			
				
## Note Model is now just the OpLines				
	for(df in 1:(length(Model)) )  {			
		if(p<Pages) {		
			nextPage<-match(p+1,Model[[df]][,1])	
			nextPageStarts<-c(nextPageStarts,nextPage)	
		}else{		
			## this makes the last page entry work	
			nextPage=length(Model[[df]][,1])+1 	
		}		
		LengthsVec<-c(LengthsVec,nextPage-thisPageStarts[df])		
		DurationsVec<-c(DurationsVec,Model[[df]]$Duration[thisPageStarts[df]:(nextPage-1)])		
		TimesVec<-c(TimesVec,Model[[df]]$Time[thisPageStarts[df]:(nextPage-1)])		
	}			
				
## this is the call to the C++ dll in the stosim library				
  fun_DF<-.Call("DetailOpLines",TimesVec, DurationsVec, LengthsVec, Names , PACKAGE="stosim")				
			
	  RcppDF<-fun_DF			
  			
	  names(RcppDF)<-c(names(RcppDF)[1:2],Names)			
				
	Page<-rep(p,length(RcppDF[,1]))			
	PageCol<-data.frame(Page)			
	RcppDF<-cbind(PageCol,RcppDF)			
				
	OutputDF<-rbind(OutputDF,RcppDF)			
				
				
	if(p==1)  { 
		if(ProgRpt==TRUE)  {	
			oneCycle<-proc.time()			
			TimeTest<-(oneCycle[3]-startTime[3])*Pages	
			  ## still need to find the correct TimeTest value				
			  if(TimeTest > .5)  {				
			  pb <- tkProgressBar(title = "MultiTrainWithInventory Progress", min = 0,				
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
				
 return(OutputDF)  				
}				
