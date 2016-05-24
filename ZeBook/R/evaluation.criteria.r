################################################################################
# "Working with dynamic models for agriculture"
# Daniel Wallach (INRA), David Makowski (INRA), James W. Jones (U.of Florida),
# Francois Brun (ACTA)
# version : 2013-03-25
################################ FUNCTIONS #####################################
#' @title Calcule multiple goodness-of-fit criteria
#' @description This function is depreciated and will be remove from the package in future versions. Please use goodness.of.fit
#' @param Ypred : prediction values from the model
#' @param Yobs : observed values
#' @param draw.plot : draw evaluation plot
#' @return data.frame with the different evaluation criteria
#' @export
#' @examples
#' # observed and simulated values
#' obs<-c(78,110,92,75,110,108,113,155,150)
#' sim<-c(126,126,126,105,105,105,147,147,147)
#' evaluation.criteria(sim,obs,draw.plot=TRUE)
evaluation.criteria=function(Ypred,Yobs,draw.plot=FALSE){
  warning("this function is depreciated, please use now goodness.of.fit")
  #Yobs is the vector of observed values
  #Ypred is the vector of predicted values
  # we keep only Ypred and Yobs where both are not NA
  select=!is.na(Ypred)&!is.na(Yobs)
  Ypred=Ypred[select]
  Yobs=Yobs[select]
  #Nobs=N  is the number of values
  Nobs<-length(Yobs)
  N<-length(Ypred)

  # Calculate and print the average of the observed values, the average of the predicted values and model bias
  mean.Yobs<-mean(Yobs,na.rm=TRUE)
  mean.Ypred<-mean(Ypred,na.rm=TRUE)
  bias<-mean.Ypred-mean.Yobs

  # Calculate and print  the standard deviation of the observed values and of the predicted values.
  # The standard deviation you need is (var(x)*(N-1)/N)^0.5 where var(x) is the R instruction that calculates sample variance.
  var.Yobs<- var(Yobs)*(Nobs-1)/Nobs
  std.Yobs<-var.Yobs^0.5
  std.Yobs
  var.Ypred<- var(Ypred)*(N-1)/N
  std.Ypred<-var.Ypred^0.5
  std.Ypred

  #Calculate and print MSE = the mean squared error and RMSE = root mean squared error.
  # unit = unit(Y)^2
  SSE<- sum((Yobs-Ypred)^2)
  SSE
  MSE<- SSE/Nobs
  MSE
  RMSE=MSE^0.5
  RMSE

  # Calculate the decomposition of MSE into 3 terms.
  # The first term is squared bias.
  bias.Squared<-bias^2
  bias.Squared
  #The second term  is SDSD ­ squared difference between the standard deviations of observed and predicted values
  SDSD<-(std.Yobs-std.Ypred)^2
  SDSD

  #The third term is the remainder term
  covBYobs<-cov(Ypred,Yobs)*(Nobs-1)/Nobs
  r<-covBYobs/(var.Yobs^0.5*var.Ypred^0.5)
  r
  LCS<-2*var.Yobs^0.5*var.Ypred^0.5*(1-r)
  LCS

  # Verify that MSE=bias.Squared+SDSD+LCS
  MSE==bias.Squared+SDSD+LCS


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


  # Calculate and print modeling efficiency
  # other possibility : EF2<-1-MSE/var.Yobs
  EF<- 1 -  sum((Yobs-Ypred)^2)/sum((Yobs-mean.Yobs)^2)
  EF
  return(data.frame(Nobs=Nobs,mean.Yobs=mean.Yobs,mean.Ypred=mean.Ypred,std.Yobs=std.Yobs,std.Ypred=std.Ypred,SSE=SSE,MSE=MSE,RMSE=RMSE,r=r,bias.Squared=bias.Squared,SDSD=SDSD,LCS=LCS,EF=EF))
}

################################################################################
# End of file
