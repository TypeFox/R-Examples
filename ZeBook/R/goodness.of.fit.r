################################################################################
# "Working with dynamic models for agriculture"
# Daniel Wallach (INRA), David Makowski (INRA), James W. Jones (U.of Florida),
# Francois Brun (ACTA)
# version : 2013-03-17
################################ FUNCTIONS #####################################
#' @title Calcule multiple goodness-of-fit criteria
#' @param Ypred : prediction values from the model
#' @param Yobs : observed values
#' @param draw.plot : draw evaluation plot
#' @return data.frame with the different evaluation criteria
#' @export
#' @examples 
#' # observed and simulated values
#' obs<-c(78,110,92,75,110,108,113,155,150)
#' sim<-c(126,126,126,105,105,105,147,147,147)
#' goodness.of.fit(obs,sim,draw.plot=TRUE)
goodness.of.fit<-function(Yobs,Ypred,draw.plot=FALSE){
  #Yobs is the vector of observed values
  #Ypred is the vector of predicted values
  Nobs<-length(Yobs)
  Npred<-length(Ypred)

  # return error if Yobs and Ypred aren't of same length
  if (Npred != Nobs)  stop("Ypred and Yobs have different lengths")

  # we keep only Ypred and Yobs where both are not NA. 
  select=!is.na(Ypred)&!is.na(Yobs)
  Ypred=Ypred[select]
  Yobs=Yobs[select]
  N<-length(Yobs)


  # Average of the observed values, the average of the predicted values and model bias
  mean.Yobs<-mean(Yobs)
  mean.Ypred<-mean(Ypred)
  bias<-mean.Yobs-mean.Ypred

  # MSE =  mean squared error  RMSE = root mean squared error.
  # RRMSE = relative root mean squared error (divide by average of observations)
  SSE<- sum((Yobs-Ypred)^2)
  MSE<- SSE/N
  RMSE=MSE^0.5
  RRMSE<-RMSE/mean.Yobs

  # MAE = mean absolute error. RMAE = relative MAE (divide by mean of observations)
  # RMAEP = realtaive MAE (divide each error by observed value)
  MAE<-mean(abs(Yobs-Ypred))
  RMAE<-MAE/mean(abs(Yobs))
  RMAEP<-mean(abs(Yobs-Ypred)/abs(Yobs))

  # modeling efficiency
  # other possibility : EF2<-1-MSE/var.Yobs
  EF<- 1 -  sum((Yobs-Ypred)^2)/sum((Yobs-mean.Yobs)^2)

  # Willmott index
  denom1<-abs(Ypred-mean(Yobs))
  denom2<-abs(Yobs-mean(Yobs))
  index<-1-sum((Yobs-Ypred)^2)/sum((denom1+denom2)^2)
 
  # Calculate the decomposition of MSE into bias², SDSD, LCS
  # The first term is squared bias.
  bias.Squared<-bias^2
  #The second term  is SDSD ­ squared difference between the standard deviations of observed and predicted values
  SDSD<-(sd(Ypred)-sd(Yobs))^2*(N-1)/N
  #The third term is the lack of correlation term
  LCS<-2*sd(Yobs)*sd(Ypred)*(1-cor(Yobs,Ypred))*(N-1)/N

  # decomposition of MSE into bias2, NU, LC
  # based on regression of observed values on simulated values
  bias2<-bias^2
  bObsPred<-cov(Yobs,Ypred)/var(Ypred)
  NU<-(1-bObsPred)^2*var(Ypred)*(N-1)/N
  LC<-(1-cor(Yobs,Ypred)^2)*var(Yobs)*(N-1)/N

  # decomposition of MSE into systematic and unsystematic parts
  # based on regression of simulated values on observed values
  bPredObs<-cov(Yobs,Ypred)/var(Yobs)
  intercept<-mean(Ypred)-bPredObs*mean(Yobs)
  newSim<-intercept+bPredObs*Yobs
  MSEs<-mean((newSim-Yobs)^2)
  MSEu<-mean((Ypred-newSim)^2)


   if( draw.plot){
    par(mfrow=c(3,1), cex=1, mar=c(4.5, 4.2, 0.5, 1.2), cex.main=0.75 )
    # observation versus prediction
    plot(Ypred,Yobs, xlim=range(c(Ypred,Yobs)),ylim=range(c(Ypred,Yobs)))
    abline(b=1,a=0)
    # residus versus prediction
    plot(Ypred,Yobs-Ypred)
    abline(a=0,b=0)
    # Plot the MSE decomposition using the instruction below
    barplot(c("MSE"=MSE,"Bias^2"=bias.Squared,"SDSD"=SDSD,"LCS"=LCS),main="decomposition of MSE", ylab="unit^2")
  }


  return(data.frame(N=N,mean.Yobs=mean.Yobs,mean.Ypred=mean.Ypred,
	bias=bias,sd.Yobs=sd(Yobs),sd.Ypred=sd(Ypred),
	MSE=MSE,RMSE=RMSE,relative.RMSE=RRMSE, 
   	MAE=MAE, relative.RMAE=RMAE, relative.MAE.prime=RMAEP, EF=EF,
	Willmott.index=index, 
	bias.squared=bias.Squared,SDSD=SDSD,LCS=LCS,
	bias.squared.again=bias.Squared, NU=NU, LC=LC,
	MSE.systematic=MSEs, MSE.unsystematic=MSEu))
}


#' @title TO COMPLETE Computation of threshold.measures
#' @param Ypred : prediction values from the model
#' @param Yobs : observed values
#' @param p : TO COMPLETE
#' @param d : TO COMPLETE
#' @param units : units
#' @return data.frame with the different evaluation criteria
#' @export
#' @examples 
#' # observed and simulated values
#' obs<-c(78,110,92,75,110,108,113,155,150)
#' sim<-c(126,126,126,105,105,105,147,147,147)
#' threshold.measures(obs,sim,80,1.0)
threshold.measures<-function(Yobs,Ypred,p,d,units=""){

  #Yobs is the vector of observed values
  #Ypred is the vector of predicted values
  Nobs<-length(Yobs)
  Npred<-length(Ypred)

  # return error if Yobs and Ypred aren't of same length
  if (Npred != Nobs)  stop("Ypred and Yobs have different lengths")

  # we keep only Ypred and Yobs where both are not NA. 
  select=!is.na(Ypred)&!is.na(Yobs)
  Ypred=Ypred[select]
  Yobs=Yobs[select]
  N<-length(Yobs)

  # find error such that at least p% of |obs-sim| are <= that value
  sortErr<-sort(abs(Yobs-Ypred))
  TDIp<-sortErr[ceiling(p/100*N)]
  print(paste("at least ",p," % of errors are below ",TDIp,units," in absolute value" ))

  # find percentage of absolute errors <=d
  CP<-100*sum(sortErr<=d)/N
  print(paste(CP," % of errors have abs value less than or equal to ",d,units))
	
  return(data.frame(TDIp=TDIp,CP=CP))
}

# end of file
