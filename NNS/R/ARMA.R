#' VN ARMA
#'
#' Autoregressive model incorporating nonlinear regressions of component series.
#'
#' @param variable Variable
#' @param h Number of periods to forecast, defaults to 1.
#' @param Training_set Sets the number of observations from 1:xx to monitor performance of forecast over in-sample range. Defaults to NULL to use all observations.
#' @param Seasonal_Factor Automatically selects the best seasonal lag from the seasonality test.  Defaults to TRUE.  To use weighted average of all seasonal lags set to FALSE.
#' @param Negative_Values If the variable can be negative, set to TRUE.  Defaults to FALSE.
#' @param Linear To use a linear regression of the component series, defaults to TRUE.  To use a nonlineaer regression, set to FALSE.
#' @param Dynamic To update the seasonal factor with each forecast point, set to TRUE.  The default is FALSE to keep the original seasonal factor from the inputted variable for all forecasts.
#' @keywords Autoressive model
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x<-rnorm(100)
#' \dontrun{VN.ARMA(x)}
#' @export



# Autoregressive Model
VN.ARMA <- function(variable,h=1,Training_set = NULL, Seasonal_Factor = TRUE ,Negative_Values = FALSE, Linear = TRUE, Dynamic = FALSE){

  original.original.variable = variable

  if(!is.null(Training_set)){
  variable = variable[1:Training_set]
  original.variable = variable[1:Training_set]
  }

  else {
    Training_set = length(variable)
    variable = variable
    original.variable = variable
    }

# Seasonality test
  output <- vector("numeric", length(variable)/4)
  instances <- vector("numeric", length(variable)/4)
  Estimates = numeric()


if(Dynamic == TRUE){
   for (j in 0:(h-1)){

    for (i in 1:floor((length(variable)/4))){

    if (abs(sd(variable[seq(length(variable),1,-i)])/mean(variable[seq(length(variable),1,-i)])) <
        abs(sd(variable)/mean(variable))){

      instances[i] <- i

      output[i]<- (abs(sd(variable[seq(length(variable),1,-i)])/mean(variable[seq(length(variable),1,-i)])))

    } else{instances[i] <- 0
    output[i]<- 0
    } }

  n<- rep(abs(sd(variable)/mean(variable)),length(instances[instances>0]))

  if(sum(instances[instances>0])==0) {lag = 1
  Weights = 1}

  if(sum(instances[instances>0])>0){
  M<- matrix(c(instances[instances>0], output[output>0],n),
             nrow=length(instances[instances>0]),  byrow= FALSE)


  Observation.sum = sum(1/M[,1])
  Observation.weighting = (1/M[,1])

  Lag.sum = sum(1/M[,2])
  Lag.weighting = (1/M[,2])

  Weights = (Lag.weighting+Observation.weighting) / (Lag.sum+Observation.sum)

  if(length(instances[instances>0])>0) {M[which.min(M[,2]),1]} else {1}

# Determine lag from seasonality test
  if(Seasonal_Factor==TRUE){lag = if(length(instances[instances>0])>0) {M[which.min(M[,2]),1]} else {1}}

  if(Seasonal_Factor==FALSE){lag<- M[,1]}
  a=length(variable)
}
  if(length(lag)==1){par(mfrow=c(2,1))} else {par(mfrow=c(1,1))}

  a=length(original.variable)
  aa=length(variable)


# Generate vectors for 1:lag
  Component.series = list()
  Component.index = list()

  for (i in 1:length(lag)){
    Component.series[[paste('Series.',i,sep="")]] <- numeric()
    Component.index[[paste('Index.',i,sep="")]] <- numeric()
    }

  for (i in 1:length(lag)){
    Component.series[[i]] = rev(variable[seq(aa+1,1,-lag[i])])
    Component.index[[i]] = (1:length(na.omit(Component.series[[i]])))
    Component.series[[i]] = na.omit(Component.series[[i]])
    }

# Regression on Component Series
  Regression.Estimates = numeric()
  Coefficients = numeric()
  reg.order = ceiling(log(aa,4))+0

    if(Linear==FALSE){
   for (i in 1:length(lag)){
    Regression.Estimates[i]=VN.ARMA.reg(Component.index[[i]],Component.series[[i]],point.est = (length(Component.series[[i]])+1),order=reg.order)
      }
    }

  if(Linear==TRUE){
    if(Negative_Values==FALSE){
  for (i in 1:length(lag)){
   Regression.Estimates[i] =  max(0,coef(lm(Component.series[[i]]~Component.index[[i]]))[1])+(coef(lm(Component.series[[i]]~Component.index[[i]]))[2]*(length(Component.series[[i]])+1))
      }
    }


      if(Negative_Values==TRUE){
        for (i in 1:length(lag)){
          Regression.Estimates[i] = (coef(lm(Component.series[[i]]~Component.index[[i]]))[1])+(coef(lm(Component.series[[i]]~Component.index[[i]]))[2]*(length(Component.series[[i]])+1))
        }
      }

  }

  Estimates[j+1] = sum(Regression.Estimates*Weights)
    variable = c(variable,(Estimates[j+1]))

  }}  #Dynamic {}

  else {


  for (i in 1:(length(original.variable)/4)){

    if (abs(sd(original.variable[seq(length(original.variable),1,-i)])/mean(original.variable[seq(length(original.variable),1,-i)])) <
        abs(sd(original.variable)/mean(original.variable))){

      instances[i] <- i

      output[i]<- (abs(sd(original.variable[seq(length(original.variable),1,-i)])/mean(original.variable[seq(length(original.variable),1,-i)])))

    }
    else{instances[i] <- 0
    output[i]<- 0
    } }



  n<- rep(abs(sd(original.variable)/mean(original.variable)),length(instances[instances>0]))
  if(sum(instances[instances>0])==0) {lag = 1
  Weights = 1}

  if(sum(instances[instances>0])>0){
  M<- matrix(c(instances[instances>0], output[output>0],n),
             nrow=length(instances[instances>0]),  byrow= FALSE)

  Observation.sum = sum(1/M[,1])
  Observation.weighting = (1/M[,1])

  Lag.sum = sum(1/M[,2])
  Lag.weighting = (1/M[,2])

  Weights = (Lag.weighting+Observation.weighting) / (Lag.sum+Observation.sum)

   {M[which.min(M[,2]),1]}

  # Determine lag from seasonality test

  if(Seasonal_Factor==TRUE){lag = M[which.min(M[,2]),1]}

  if(Seasonal_Factor==FALSE){lag<- M[,1]}
  a=length(variable)
  }

  if(length(lag)==1){par(mfrow=c(2,1))} else {par(mfrow=c(1,1))}



 #if(is.null(Dynamic)) {
  for (j in 0:(h-1)){

    a=length(original.variable)
    aa=length(variable)


    # Generate vectors for 1:lag
    Component.series = list()
    Component.index = list()

    for (i in 1:length(lag)){
      Component.series[[paste('Series.',i,sep="")]] <- numeric()
      Component.index[[paste('Index.',i,sep="")]] <- numeric()
    }

    for (i in 1:length(lag)){
      Component.series[[i]] = rev(variable[seq(aa+1,1,-lag[i])])
      Component.index[[i]] = (1:length(na.omit(Component.series[[i]])))
      Component.series[[i]] = na.omit(Component.series[[i]])
    }



    # Regression on Component Series
    Regression.Estimates = numeric()
    Coefficients = numeric()
    reg.order = ceiling(log(aa,4))+0

    if(Linear==FALSE){
      for (i in 1:length(lag)){
        Regression.Estimates[i]=VN.ARMA.reg(Component.index[[i]],Component.series[[i]],point.est = (length(Component.series[[i]])+1),order=reg.order)
        }
      }

    if(Linear==TRUE){
      if(Negative_Values==FALSE){
        for (i in 1:length(lag)){
          Regression.Estimates[i] =  max(0,coef(lm(Component.series[[i]]~Component.index[[i]]))[1])+(coef(lm(Component.series[[i]]~Component.index[[i]]))[2]*(length(Component.series[[i]])+1))
        }
      }


      if(Negative_Values==TRUE){
        for (i in 1:length(lag)){
          Regression.Estimates[i] = (coef(lm(Component.series[[i]]~Component.index[[i]]))[1])+(coef(lm(Component.series[[i]]~Component.index[[i]]))[2]*(length(Component.series[[i]])+1))
        }
      }

    }

    Estimates[j+1] = sum(Regression.Estimates*Weights)
    variable = c(variable,(Estimates[j+1]))


  }}  # ELSE {}

#### PLOTTING
  par(mfrow=c(2,1))
  if(sum(instances[instances>0])>0){
  plot(instances[instances>0],output[output>0],
       xlab="Period", ylab="Coefficient of Variance", main = "Seasonality Test",
       ylim = c(0,2*abs(sd(original.variable)/mean(original.variable))),
       col=ifelse(output[output>0]==min(output[output>0]), "red", "black"),
       pch =ifelse(output[output>0]==min(output[output>0]), 19, 1))

  abline(h=abs(sd(original.variable)/mean(original.variable)), col="red",lty=5)
  text(length(instances[instances>0])/2,abs(sd(original.variable)/mean(original.variable)),adj=c(0,-0.25),
      "Variable Coefficient of Variance",col='red')
  }

  plot( original.original.variable, type = 'l',main = "Forecast",col='steelblue',  xlim=c(1,max((a+h),length(original.original.variable))), ylab="Variable", ylim=c(min(Estimates, original.original.variable),max( original.original.variable,Estimates)))

  lines((Training_set+1):(Training_set+h),Estimates,type = 'l',lty=3,col='red')
  segments(Training_set,original.variable[Training_set],Training_set+1,Estimates[1],lty=3,col='red')

  legend('topleft',bty='n',legend = c("Original", paste("Forecast ",h," period(s)",sep = "")),lty = c(1,2),col=c('steelblue','red'))
  points(Training_set,variable[Training_set],col="green",pch=18)
  points(Training_set+h,sum(Regression.Estimates*Weights),col="green",pch=18)

#  if(Negative_Values==TRUE){
 # for (i in 1:length(lag)){
 #   segments(0,coef(lm(Component.series[[i]]~Component.index[[i]]))[1],
#             length(original.variable)+h,
#           Regression.Estimates[i]
#           ,col=adjustcolor("red", alpha.f = .125), lty = 3)
#  }
#  }

#  if(Negative_Values==FALSE){
#    for (i in 1:length(lag)){
#      segments(0,max(0,coef(lm(Component.series[[i]]~Component.index[[i]]))[1]),
 #              length(original.variable)+h,
 #              max(0,Regression.Estimates[i])
#               ,col=adjustcolor("red", alpha.f = 0.125), lty = 3)
#
 #   }
#  }

return(Estimates)

}


