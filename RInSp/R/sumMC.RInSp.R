sumMC.RInSp = function (dataset) 
{
  #
  # The function  provides summary statistics on the Monte Carlo resampling simulated index distributions.
  #
  # Author: Nicola ZACCARELLI, Giorgio MANCINELLI, Dan BOLNICK
  # E-mail: nicola.zaccarelli@gmail.com,
  #         giorgio.mancinelli@unisalento.it
  #         danbolnick@mail.texas.edu
  #
  # Version: 1.0
  # Date: 10/11/2012
  #
    if (class(dataset) != "RInSp") 
        stop("The input must be an object of class RInSp from a Monte Carlo simulation")
    if (dataset$parameter == 0) stop("There are not Monte Carlo data to plot.")
   # Let's build the title for the graph and pick right x label
   if (dataset$parameter == 1) { 
       TitleGraph = expression(paste("Distribution of ", PS[i]))
       LableX = expression(PS[i])
     }
   if (dataset$parameter == 2) { 
       TitleGraph = "Distribution of E"
       LableX = "E"
     }
   if (dataset$parameter == 4) { 
       TitleGraph = "Distribution of WIC/TNW"
       LableX = "WIC/TNW"
     }
    if (dataset$parameter == 5) { 
       TitleGraph = "Distribution of WIC/TNW"
       LableX = "WIC/TNW"
     }
    hist(dataset$montecarlo[, dataset$parameter], freq = FALSE, main = TitleGraph, xlab = LableX,
         breaks= 2 * nclass.Sturges(dataset$montecarlo[, dataset$parameter]), col="grey70")
    legend(x="topright", legend=c("sample value", "95% conf. limits"), col=c("blue", "red"), lty=c(1, 2))
  # Plotting calculated value's bar
    datahist =  hist(dataset$montecarlo[, dataset$parameter], breaks= 2 * nclass.Sturges(dataset$montecarlo[, dataset$parameter]), plot=FALSE)
    segments(x0=dataset$montecarlo[1, dataset$parameter], y0=0, x1=dataset$montecarlo[1, dataset$parameter],
             y1= max(datahist$density)/2, col="blue", lty = 1, lwd = 2)
   # plot 2.5% and 97.5 bootstrapped percentile
    percboot = quantile(dataset$montecarlo[ , dataset$parameter], probs=c(0.025, 0.975), na.rm = TRUE)
    segments(x0=percboot[1], y0=0, x1=percboot[1], y1= max(datahist$density)/2, col="red", lty = 2, lwd = 2)
    segments(x0=percboot[2], y0=0, x1=percboot[2], y1= max(datahist$density)/2, col="red", lty = 2, lwd = 2)
return(summary(dataset$montecarlo[, dataset$parameter]))
}

