#'Niche overlap null models 
#'@description Create a null model for niche overlap; choices of algorithm and metric are constrained to be valid for niche null models.
#'@param speciesData a data frame in which each row is a species, each column is a resource utilization category, and the entries represent the quantity of the resource used by each species. Examples might be the amount of time a species spends foraging in different microhabitats, the biomass of different prey types, or counts of the number of times an adult female oviposits eggs on different species of a host plant.
#'@param algo the algorithm to use, must be "ra1", "ra2", "ra3", "ra4"; default is "ra3".
#'@param metric the metric used to calculate the null model: choices are "pianka", "czekanowski", "pianka_var", "czekanowski_var", "pianka_skew", "czekanowski_skew"; default is "pianka".
##'@param nReps the number of replicate null assemblages to create; default is 1000 replicates.
#'@param saveSeed TRUE or FALSE. If TRUE the current seed is saved so the simulation can be repeated; default is FALSE.
#'@param algoOpts a list containing all the options for the specific algorithm you want to use.  Must match the algorithm given in the `algo` argument.
#'@param metricOpts a list containing all the options for the specific metric you want to use.  Must match the metric given in the `metric` argument.
#'@param suppressProg TRUE or FALSE. If true, display of the progress bar in the console is suppressed; default is FALSE. This setting is useful for creating markdown documents with `knitr`.
#'@examples \dontrun{
#' ## Load MacAruthur warbler data
#' data(dataMacWarb)
#' 
#' ## Run the null model
#' warbMod <- niche_null_model(dataMacWarb,nReps=1000)
#' ## Summary and plot info
#' summary(warbMod)
#' plot(warbMod)
#' plot(warbMod,type="niche")
#'}
#'
#'@export

niche_null_model <- function(speciesData, algo = "ra3", metric = "pianka", nReps = 1000,saveSeed=FALSE,algoOpts = list(),metricOpts = list(),suppressProg = FALSE){
  
  aChoice <- c("ra1","ra2","ra3","ra4")
  mChoice<- c("pianka", "czekanowski", "pianka_var", "czekanowski_var", "pianka_skew", "czekanowski_skew")
  algo <- match.arg(algo,choices = aChoice)
  metric <- match.arg(metric,choices = mChoice)
  params <- list(speciesData = speciesData, algo = algo, metric = metric, nReps = nReps, saveSeed = saveSeed,algoOpts = algoOpts,metricOpts = metricOpts, suppressProg = suppressProg)
  output <- do.call(null_model_engine,params)
  class(output) <- "nichenullmod"
  return(output)

}


#' Generic function for calculating null model summary statistics
#' @description Takes as input a list of Null.Model.Out, with Obs, Sim, Elapsed Time, and Time Stamp values
#' @param object the null model object to print a summary of.
#' @param ... Extra parameters for summary .
#' @details The summary output includes a timestamp and complete statistics on the simulated values of the metric, including the mean, variance, and one and two-tailed 95% confidence intervals. The lower and upper tails for the observed value are given, as is the standardized effect size (SES), which is calculated as the (Metric(obs) - average(Metric(sim)))/(standard deviation(Metric(sim))). Large positive values (or negative) values indicate that the observed metric is significantly larger (or smaller) than predicted by the null model. If the distribution of errors is approximately normal, then non-significant values will fall roughly within +- two SES values. 
#' @export

summary.nichenullmod <- function(object,...)
{ 
 
  nullmodObj <- object 
  
  cat("Time Stamp: " , nullmodObj$Time.Stamp,   "\n")  
  cat("Reproducible: ",nullmodObj$Reproducible,  "\n")
  cat("Number of Replications: ",nullmodObj$nReps,  "\n")
  cat("Elapsed Time: ", nullmodObj$Elapsed.Time, "\n")
  cat("Metric: ", nullmodObj$Metric,  "\n")
  cat("Algorithm: ", nullmodObj$Algorithm,  "\n") 
  
  cat("Observed Index: ", format(nullmodObj$Obs,digits=5),  "\n")
  cat("Mean Of Simulated Index: ",format(mean(nullmodObj$Sim),digits=5),  "\n")
  cat("Variance Of Simulated Index: ",format(var(nullmodObj$Sim),digits=5),  "\n")
  cat("Lower 95% (1-tail): ",format(quantile(nullmodObj$Sim,0.05),digits=5),  "\n")
  cat("Upper 95% (1-tail): ",format(quantile(nullmodObj$Sim,0.95),digits=5), "\n")
  cat("Lower 95% (2-tail): ",format(quantile(nullmodObj$Sim,0.025),digits=5), "\n")
  cat("Upper 95% (2-tail): ",format(quantile(nullmodObj$Sim,0.975),digits=5),  "\n")
  
  #P-values
  if (nullmodObj$Obs > max(nullmodObj$Sim)) {
    cat("Lower-tail P > ",(length(nullmodObj$Sim) - 1)/length(nullmodObj$Sim),  "\n")
    cat("Upper-tail P < ",1/length(nullmodObj$Sim),  "\n")
  } else if(nullmodObj$Obs < min(nullmodObj$Sim)) {
    cat("Lower-tail P > ", 1/length(nullmodObj$Sim), "\n")
    cat("Upper-tail P < ",(length(nullmodObj$Sim) - 1)/length(nullmodObj$Sim), "\n")
  } else {
    cat("Lower-tail P = ", format(sum(nullmodObj$Obs >= nullmodObj$Sim)/length(nullmodObj$Sim),digits=5),  "\n")
    cat("Upper-tail P = ", format(sum(nullmodObj$Obs <= nullmodObj$Sim)/length(nullmodObj$Sim),digits=5), "\n")
  }

  cat(paste("Observed metric > ",sum(nullmodObj$Obs > nullmodObj$Sim)," simulated metrics",sep="") , "\n")
  cat(paste("Observed metric < ",sum(nullmodObj$Obs < nullmodObj$Sim)," simulated metrics",sep="") , "\n")
  cat(paste("Observed metric = ",sum(nullmodObj$Obs == nullmodObj$Sim)," simulated metrics",sep=""),  "\n")
  cat("Standardized Effect Size (SES): ", format((nullmodObj$Obs - mean(nullmodObj$Sim))/sd(nullmodObj$Sim),digits=5), "\n")
  
  #if(!is.null(Output.File)) close(outfile)
}




#' Niche Null Model Plot function
#' @description Plot niche overlap null model object.
#' @param x the null model object to plot.
#' @param type the type of null model plot to display.  See details for more information.
#' @param ... Other variables to be passed on to base plotting.
#' @details the valid types for the Niche Overlap module are "hist" to display a histogram of the simulated metric values, and "niche" to display the observed data matrix and one simulated matrix.
#' 
#' The "hist" plot type is common to all EcoSimR modules. The blue histogram represents the NRep values of the metric for the simulated assemblages. The red vertical line represents the metric value for the real assemblage. The two pairs of vertical dashed black lines represent the one-tailed (long dash) and two-tailed (short dash) 95% confidence exact confidence intervals of the simulated data.
#' 
#' The "niche" plot type illustrates the utilization data (observed = red, simulated = blue). Each row in the figure is a species, and each column is a resource utilization category. The area of each circle depicted is proportional to the utilization of a resoruce category by a species. Empty positions indicate a resource utilization of 0.0. The rows and columns are illustrated with the same ordering as the original data matrix.
#' 
#' @export



plot.nichenullmod <- function(x, type = "hist",...)
{
  nullmodObj <- x
  if(type == "hist"){

  opar <- par(no.readonly=TRUE)
  par(cex=1, cex.axis = 1.5,
      cex.main=1,cex.lab=1.6)
  par (mar=c(5,6,4,2)+0.1,mfrow=c(1,1))
  #------------------------------------------------------
  hist(nullmodObj$Sim, breaks=20, col="royalblue3",
       
       xlab="Simulated Metric",ylab="Frequency",main="",
       xlim=range(c(nullmodObj$Sim,nullmodObj$Obs)))
  abline(v=nullmodObj$Obs,col="red",lty="solid",lwd=2.5)
  abline(v=quantile(nullmodObj$Sim,c(0.05,0.95)),col="black",lty="dashed",lwd=2.5)
  abline(v=quantile(nullmodObj$Sim,c(0.025,0.975)),col="black",lty="dotted",lwd=2.5)
  mtext(as.character(date()),side=3,adj=1,line=3)
  }
  
  if(type == "niche"){
    opar<- par(no.readonly=TRUE)
    par(mfrow=c(2,1))
    Data <- nullmodObj$Data/rowSums(nullmodObj$Data)
    plot(rep(1:ncol(Data),times = nrow(Data)),
         rep(1:nrow(Data),each=ncol(Data)),
         xlab="Resource Category",ylab="Species",cex=10*sqrt(t(Data)/pi),col="red3",lwd=2,
         main="Observed Utilization Matrix",col.main="red3",cex.main=1.5)
  mtext(as.character(nullmodObj$Time.Stamp),side=3,adj=1,line=3)
    
    One.Null.Matrix <- nullmodObj$Randomized.Data
    One.Null.Matrix <- One.Null.Matrix/rowSums(One.Null.Matrix)
    plot(rep(1:ncol(One.Null.Matrix),times = nrow(One.Null.Matrix)),
         rep(1:nrow(One.Null.Matrix),each=ncol(One.Null.Matrix)),
         xlab="Resource Category",ylab="Species",cex=10*sqrt(t(One.Null.Matrix)/pi),col="royalblue3",lwd=2,
         main="Simulated Utilization Matrix",col.main="royalblue3",cex.main=1.5)
    par(opar)
  }
  
}
