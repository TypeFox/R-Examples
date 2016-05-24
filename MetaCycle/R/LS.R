### Authors of original R code of Lomb-Scargle: Earl F. Glynn, Jie Chen and Arcady R. Mushegian
### Email: chenj@umkc.edu
### Associated literature: Earl F. Glynn, Jie Chen and Arcady R. Mushegian. Bioinformatics. 22(3):310-6 (2006).
### Website: http://research.stowers-institute.org/mcm/efg/2005/LombScargle/R/index.htm
### Lomb-Scargle Normalized Periodogram:  
    ## "Fast Algorithm for Spectral Analysis of Unevenly Sampled Data". William H. Press and George B. Rybicki. Astrophysical Journal, 338:277-280, March 1989.
    ## Also appeared in Section 13.8, "Spectral Analysis of Unevenly Sampled Data" in Numerical Recipes in C (2nd Ed). William H. Press, et al, Cambridge University Press, 1992.
#============================================================================================      
NHorneBaliunas <- function(N)                                                 #'NHorneBaliunas' function is from Lomb-Scargle.R                           
{
  # Estimate of number of independent frequencies in Lomb-Scargle analysis based on sample size. 
  # Use this regression equation for estimating the number of independent frequencies for calls to ComputeLombScargle or ComputeAndPlotLombScargle.
  # From Horne and Baliunas, "A Prescription for Period Analysis of Unevenly Sampled Time Series",The Astrophysical Journal, 392: 757-763, 1986.
  #----------------------
  Nindependent <- trunc(-6.362 + 1.193*N + 0.00098*N^2)                       ##if N=6, -6.362 + 1.193*N + 0.00098*N^2=0.83128; trunc(0.83128)=0;
  if (Nindependent < 1)                                                       
  { Nindependent <- 1   }                                                     # Kludge for now
  return(Nindependent)
}
#-------------------------------------------------------------------------------------------
ComputeLombScargle <- function(t, h, TestFrequencies, Nindependent)           ##'ComputeLombScargle' function is from Lomb-Scargle.R 
{
  # "h" is vector of expression values for time points "t".
  # SpectralPowerDensity will be evaluated at given TestFrequencies.
  # Nindepedent of the TestFrequencies are assumed to be independent.
  #----------------------
  stopifnot( length(t) == length(h) )                                        
  if (length(t) > 0)
  {
    Nyquist <- 1 / (2 * ( (max(t) - min(t) )/ length(t) ) )                   ##get the value of Nyquist(the highest frequency for which the unevenly spaced data may be evaluated)
    hResidual    <- h - mean(h)                                               
    SpectralPowerDensity <- rep(0, length(TestFrequencies))       
    for (i in 1:length(TestFrequencies))
    {
      ##The values (eg. Omega, TwoOmegaT, Tau, OmegaTMinusTau, and SpectralPowerDensity) are calculated using the formula mentioned in the Lomb-Scargle paper.
   	  Omega       <- 2*pi*TestFrequencies[i]                                  
      TwoOmegaT   <- 2*Omega*t          
      Tau         <- atan2( sum(sin(TwoOmegaT)) , sum(cos(TwoOmegaT)) ) / (2*Omega)            ##for positive arguments atan2(y, x) == atan(y/x)
      OmegaTMinusTau <- Omega * (t - Tau)     
      SpectralPowerDensity[i] <- (sum (hResidual * cos(OmegaTMinusTau))^2) / sum( cos(OmegaTMinusTau)^2 ) + 
                                (sum (hResidual * sin(OmegaTMinusTau))^2) / sum( sin(OmegaTMinusTau)^2 )
    }
    # The "normalized" spectral density refers to the variance term in the denominator. 
	# With this term the SpectralPowerDensity has an exponential probability distribution with unit mean.
    SpectralPowerDensity <- SpectralPowerDensity / ( 2 * var(h) )             
    Probability <- 1 - (1-exp(-SpectralPowerDensity))^Nindependent            ##get the probability
    PeakIndex    <- match(max(SpectralPowerDensity), SpectralPowerDensity)    ##get the the peak index
    # Note:  Might merit more investigation when PeakIndex is the first point.
    PeakPeriod <- 1 / TestFrequencies[PeakIndex]                              ##get the period length corresponding to the peak index with largest normalized spectral power density value
    PeakPvalue <- Probability[PeakIndex]                                      ##get the corresponding p-value to the peak index with largest normalized spectral power density value
  } else {                                                                    
    # Time series has 0 points
    Nyquist     <- NA
    Probability <- 1.0
    PeakIndex   <- NA
    SpectralPowerDensity <- NA
    PeakPeriod  <- NA
    PeakPvalue  <- 1.0
  }
  return( list( t=t,h=h,Frequency=TestFrequencies,Nyquist=Nyquist,
                SpectralPowerDensity=SpectralPowerDensity,Probability=Probability,
                PeakIndex=PeakIndex,PeakSPD=SpectralPowerDensity[PeakIndex],
				PeakPeriod=PeakPeriod,PeakPvalue=PeakPvalue,N=length(h),Nindependent=Nindependent) )
}
#-------------------------------------------------------------------------------------------
ComputeAndPlotLombScargle <- function(t, h, TestFrequencies, Nindependent)    ##'ComputeAndPlotLombScargle' function is from Lomb-Scargle.R
{
  h <- h[ order(t) ]                                                          # order expression values in ascending time order   
  t <- t[ order(t) ]                                                          # order time points in ascending time order
  t2 <- t[!(is.na(h) | is.nan(h))]                                            # Remove missing (na) and not-a-number (NaN) expression values    
  h2 <- h[!(is.na(h) | is.nan(h))]         
  h2 <- h2[ order(t2) ]                                                       # Order by time                                                     
  t2 <- t2[ order(t2) ]                                                       
  N <- length(h2)            
  LS <- ComputeLombScargle(t2, h2, TestFrequencies, Nindependent)  
  LS$t.all     <- t                                                           
  LS$h.all     <- h                                                           
  if (N > 5)                                            
  {
    # Compute loess smoothed curve and find peak (assume only one for now)
    loess.fit <- loess(h ~ t, data.frame(t=t, h=h))                        
    h.loess   <- predict(loess.fit, data.frame(t=t))                        
    h.peak    <- optimize(function(t, model)  predict(model, data.frame(t=t)), c(min(t),max(t)),maximum=TRUE,model=loess.fit)
						                                                      ##get the maximum value of smoothed profile value and its corresponding time points value
    LS$h.loess   <- h.loess 
    LS$h.peak    <- h.peak 
  } else {                                                                           
    LS$h.loess <- NULL     
    LS$h.peak  <- list(maximum=NaN,objective=NaN)
  }
  return(LS)
}
###======================================================================================================================================
runLS <- function(indata,LStime,minper=20,maxper=28,releaseNote=TRUE)
{
	if (releaseNote)  {
		cat("The LS is in process from ", format(Sys.time(), "%X %m-%d-%Y"),"\n");
	}
	RawData <- indata;
	outID <- as.character(RawData[,1]);
	Expression <- data.matrix(RawData[,2:ncol(RawData)]);
	Time <- LStime;
	stopifnot( length(outID)   == nrow(Expression) );                               
	#-----------------------
	N <- ncol(Expression);
	stopifnot(length(Time) == N);
	#for keeping the 'TestFrequencies' same when adding columns of  missing values to the dataset
	emptyNum <- 0;
	for (i in 1:ncol(Expression))
	{
		if (all(is.na(Expression[,i])))
		{	emptyNum <- emptyNum + 1;  }
	}
	M <- 4*(N - emptyNum);
	MinFrequency <- 1/maxper;
	MaxFrequency <- 1/minper;
	if ( (MaxFrequency > 1/(2*mean(diff(Time)))) & (releaseNote) )
	{  cat("MaxFrequency may be above Nyquist limit.\n"); }
	TestFrequencies <- MinFrequency +
					  (MaxFrequency - MinFrequency) * (0:(M-1) / (M-1))
	#-----------------------
	header<-c("CycID","PhaseShift","PhaseShiftHeight","PeakIndex",
	          "PeakSPD","Period","p","N","Nindependent","Nyquist");
	LSoutM<-header;
	for (j in 1:length(outID))                                                     
	{
		if (.Platform$OS.type == "windows")                                       
		{ flush.console(); }                                                       # display immediately in Windows 
		Nindependent <- NHorneBaliunas(sum(!is.na(Expression[j,]))); 
		LS <- ComputeAndPlotLombScargle(Time, Expression[j,], TestFrequencies, Nindependent);
		LSoutM <- rbind(LSoutM,c(outID[j], LS$h.peak$maximum,LS$h.peak$objective, LS$PeakIndex,LS$PeakSPD, 
							  LS$PeakPeriod, LS$PeakPvalue, LS$N, LS$Nindependent,  LS$Nyquist));
	}
	LSoutM<-LSoutM[2:nrow(LSoutM),];
    colnames(LSoutM) <- header;
	bhq <- p.adjust(as.numeric(LSoutM[,"p"]),"BH");
	LSoutM <- cbind(LSoutM,bhq);	       
	dimnames(LSoutM)<-list("r"=1:length(outID),"c"=c(header,"BH.Q"));
	if (releaseNote)  {
		cat("The analysis by LS is finished at ",format(Sys.time(), "%X %m-%d-%Y"),"\n");
	}
	return(LSoutM);
}
#============================================================================================  
