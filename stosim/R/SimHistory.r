   SimHistory<-function(Model,SimulationYears = 2000,SimulationYearsPerPage = 1000,
   ProgRpt=FALSE) {      
  
  ShowProgress=FALSE
  thisRNGkind="Marsaglia-Multicarry"
				
	## validate input dataframe for Type1 Model			
	Type1_Names<-c("OpLine","EventID","FD","FP1","FP2","FP3","RD","RP1","RP2","RP3","Seed")			
	Input_Names<-names(Model)			
				
	if(length(Input_Names)!=length(Type1_Names))  {			
	stop("Error in stosim:    incorrect input dataframe length to SimHistory")			
	} else {			
	for(n in 1:length(Type1_Names)) {			
	if(Input_Names[n]!=Type1_Names[n]) stop("Error in stosim:   incorrect input dataframe names to SimHistory") }			
	}			
				
	if(length(Model[1])>1) {			
	for(n in 1:(length(Model[,1])-1))  {			
	if( Model[n,1] != Model[n+1,1] )  stop("Error in stosim:   This is not a Type1 model")			
	} }			
				
	Pages=as.integer(SimulationYears/SimulationYearsPerPage)			
	if(Pages<1) {			
	stop("Error in stosim:  SimulationYears are less than SimulationYearsPerPage")			
	}			

  if(SimulationYearsPerPage>2000)  {
    stop("Error in stosim:  Accuracy lost over 2000 years per page")
  }
  
	## Extract vectors to pass to Rcpp			
	OpLine_Vec<-Model[,1]			
	Event_ID_Vec<-Model[,2]			
	FD_DF<-NULL			
	## Parse FD to convert type to a code			
				
	for(x in 1:length(Model[,3]))  {			
				
	if (Model[x,3]=="E")  FD_DF<-rbind(FD_DF,1)			
	if (Model[x,3]=="N")  FD_DF<-rbind(FD_DF,2)			
	if (Model[x,3]=="W")  FD_DF<-rbind(FD_DF,3)			
	## if (Model[x,3]=="L")  FD_DF<-rbind(FD_DF,4)			
	if(length(FD_DF[,1])<x) stop("Error in stosim:   FD parameter not recognized")			
	}			
	FD_Vec<-as.vector(FD_DF)			
	FP1_Vec<-Model[,4]			
	FP2_Vec<-Model[,5]			
	FP3_Vec<-Model[,6]			
	RD_DF<-NULL			
	## Parse RD to convert type to a code			
				
	for(x in 1:length(Model[,7]))  {			
				
	##if (Model[x,7]=="E")  RD_DF<-rbind(RD_DF,1)			
	if (Model[x,7]=="N")  RD_DF<-rbind(RD_DF,2)			
	if (Model[x,7]=="W")  RD_DF<-rbind(RD_DF,3)			
	if (Model[x,7]=="L")  RD_DF<-rbind(RD_DF,4)			
	if(length(RD_DF[,1])<x) stop("Error in stosim:   RD parameter not recognized")			
	}			
	RD_Vec<-as.vector(RD_DF)			
	RP1_Vec<-Model[,8]			
	RP2_Vec<-Model[,9]			
	RP3_Vec<-Model[,10]			
	Seed_Vec<-as.integer(Model[,11])			
				
	OutputDF=NULL			
## even proc.time has a problem on CRAN example run
if(ProgRpt==TRUE)  {			
	startTime<-proc.time()	
}	  
  
## the main loop through pages				
for(p in 1:Pages)  {				
				
	if(p>1) { Seed_Vec<-Seed_Vec +3 }			
				
				
	## use this set.seed to establish the ring.kind for the Rcpp call			
	RNGkind(kind=thisRNGkind)			
				
  ## this is the call to the C++ code in the stosim library				
  fun_out<-.Call("SimulationHistory",				
    OpLine_Vec, Event_ID_Vec, FD_Vec, FP1_Vec, FP2_Vec,				
    FP3_Vec, RD_Vec, RP1_Vec, RP2_Vec, RP3_Vec, Seed_Vec,				
    SimulationYearsPerPage, PACKAGE="stosim")				
				
	Page<-rep(p,length(fun_out[,1]))			
	PageCol<-data.frame(Page)			
	fun_out<-cbind(PageCol,fun_out)			
				
	OutputDF<-rbind(OutputDF,fun_out)	
  
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
	setTkProgressBar(pb, p, title=paste( round(p/Pages*100, 0),"% done"))			
	}
  
## return to main loop through pages  
}
   
   
	if(ShowProgress==TRUE)  {			
	close(pb)			
	}			
    return(OutputDF)				
    }	