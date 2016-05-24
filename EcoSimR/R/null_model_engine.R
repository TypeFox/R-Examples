#'Run null model
#'@description This function drives all the different kinds of null models that can be run. It is the underlying engine.
#'@param speciesData a dataframe for analysis that is compatable with the metrics and algorithms used. 
#'@param algo the algorithm used to randomize the data.
#'@param metric the metric used to quantify pattern in the data.
#'@param nReps the number of null assemblages to simulate.
#'@param saveSeed Save the existing random seed to allow the user to reproduce the exact model results. The default value is FALSE, in which case the random number seed that is created is not stored in the output. 
#'@param algoOpts a list containing options for a supplied alogrithm.
#'@param metricOpts a list containing options for a supplied metric.
#'@param type The type of null model being run. If the null model is intended to be used with one of the existing modules, the type should be "size","niche", or "cooc". If the user is creating an entirely new null model, type should be set to NULL (the default value).
#'@param suppressProg TRUE or FALSE. If true, display of the progress bar in the console is suppressed; default is FALSE. This setting is useful for creating markdown documents with `knitr`.
#'@examples \dontrun{
#' # User defined function
#' 
#'
#'}
#'@export


null_model_engine <- function(speciesData, algo, metric, nReps = 1000, saveSeed = FALSE, algoOpts = list(), metricOpts = list(),type=NULL,suppressProg=FALSE)
{
  if(suppressProg){
    pb <- txtProgressBar(min = 0, max = nReps, style = 3, file = stderr())
  } else{
  pb <- txtProgressBar(min = 0, max = nReps, style = 3)
  }
  ## Set the seed
  if(saveSeed){
    randomSeed <- .Random.seed
  } else {
    randomSeed <- NULL
  }
 
  ## Convert to matrix for type consistency
  if(!is.matrix(speciesData)){ speciesData <- as.matrix(speciesData)}
  
  ### Check for row names hidden in the data frame and automagically strip them.
  
  if(suppressWarnings(is.na(as.numeric(speciesData[2,1])))){
    speciesData <- speciesData[,-1] 
    class(speciesData) <- "numeric"
  }

  algoF <- get(algo)
  metricF <- get(metric)
  
  ### Error check for input functions
  if(!grepl("speciesData",names(formals(algoF)[1]))){
    stop("Please enter a valid algorithm with 'speciesData' as the first parameter")
  }
  if(!grepl("m",names(formals(metricF)[1])) && nchar(names(formals(metricF)[1])) == 1 ){
    stop("Please enter a valid metric with 'm' as the first parameter")
  }
    
  
  
  startTime <- Sys.time()
  
  
  if (nReps < 2) nReps <- 2
  
  sim <- rep(NA,nReps)

  algoOpts[["speciesData"]] <- speciesData
  metricOpts[["m"]] <- speciesData
  obs <- as.vector(do.call(metricF,metricOpts))
  
  
  for(i in 1:nReps){
    m <- do.call(algoF,algoOpts)
    metricOpts[["m"]] <- m
    sim[i] <- do.call(metricF,metricOpts)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  ## Save final matrix for plotting
  finalRandomData <- m
  
  endTime <- Sys.time()
  elapsedTime <- format(endTime-startTime,digits=2)
  timeStamp <- date()
  
  ### Reverse engineers the naming for consistent output
  ### Be sure to update the code below if new algos and metrics are added
  

  
  nullModelOut <- list(Obs=obs,Sim=sim, Elapsed.Time=elapsedTime, Time.Stamp=timeStamp,Metric = metric, Algorithm = algo, nReps = nReps, 
                       Reproducible = saveSeed,RandomSeed = randomSeed, Data = speciesData,Randomized.Data = finalRandomData)
  
  if(is.null(type)){
  class(nullModelOut) <- "nullmod"
} else if (type %in% c("niche","cooc","size")){
  class(nullModelOut) <- paste(type,"nullmod",sep="")
  
  
}
  
  return(nullModelOut)
  
}



#' Generic function for calculating null model summary statistics
#' @description Takes as input a list of Null.Model.Out, with Obs, Sim, Elapsed Time, and Time Stamp values
#' @param object the null model object to print a summary of.
#' @param ... Extra parameters for summary.
#' @details The summary output includes a timestamp and complete statistics on the simulated values of the metric, including the mean, variance, and one and two-tailed 95% confidence intervals. The lower and upper tails for the observed value are given, as is the standardized effect size (SES), which is calculated as the (Metric(obs) - average(Metric(sim)))/(standard deviation(Metric(sim))). Large positive values (or negative) values indicate that the observed metric is significantly larger (or smaller) than predicted by the null model. If the distribution of errors is approximately normal, then non-significant values will fall roughly within +- two SES values. 
#' @export

summary.nullmod <- function(object,...)
{ 
  nullmodObj <- object
  #if (!is.null(Output.File)) outfile <- file(p$Output.File, "w") else outfile <-""
  
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
    cat("Lower-tail P < ",(length(nullmodObj$Sim) - 1)/length(nullmodObj$Sim),  "\n")
    cat("Upper-tail P > ",1/length(nullmodObj$Sim),  "\n")
  } else if(nullmodObj$Obs < min(nullmodObj$Sim)) {
    cat("Lower-tail P > ", 1/length(nullmodObj$Sim), "\n")
    cat("Upper-tail P < ",(length(nullmodObj$Sim) - 1)/length(nullmodObj$Sim), "\n")
  } else {
    cat("Lower-tail P = ", format(sum(nullmodObj$Obs >= nullmodObj$Sim)/length(nullmodObj$Sim),digits=5),  "\n")
    cat("Upper-tail P = ", format(sum(nullmodObj$Obs <= nullmodObj$Sim)/length(nullmodObj$Sim),digits=5), "\n")
  }
  
  cat(paste("Observed metric > ",sum(nullmodObj$Obs > nullmodObj$Sim)," simulated metrics",sep="") , "\n")
  cat(paste("Observed metric < ",sum(nullmodObj$Obs < nullmodObj$Sim)," simulated metrics",sep="")  ,"\n")
  cat(paste("Observed metric = ",sum(nullmodObj$Obs == nullmodObj$Sim)," simulated metrics",sep="") , "\n")
  cat("Standardized Effect Size (SES): ", format((nullmodObj$Obs - mean(nullmodObj$Sim))/sd(nullmodObj$Sim),digits=5), "\n")
  
  #if(!is.null(Output.File)) close(outfile)
}


#' plot a histogram null model
#' @description Plot a null model object.
#' @param x the null model object to plot.
#' @param ... Other variables to be passed on to base plotting.
#' @details The "hist" plot type is common to all EcoSimR modules. The blue histogram represents the NRep values of the metric for the simulated assemblages. The red vertical line represents the metric value for the real assemblage. The two pairs of vertical dashed black lines represent the one-tailed (long dash) and two-tailed (short dash) 95% confidence exact confidence intervals of the simulated data.
#' 
#' @export



plot.nullmod <- function(x,...)
{
  nullmodObj <- x
    par(mfrow=c(1,1))
    opar <- par(no.readonly=TRUE)
    par(cex=1, cex.axis = 1.5,
        cex.main=1,cex.lab=1.6)
    par (mar=c(5,6,4,2)+0.1)
    hist(nullmodObj$Sim, breaks=20, col="royalblue3",
         
         xlab="Simulated Metric",ylab="Frequency",main="",
         xlim=range(c(nullmodObj$Sim,nullmodObj$Obs)))
    
    abline(v=nullmodObj$Obs,col="red",lty="solid",lwd=2.5)
    abline(v=quantile(nullmodObj$Sim,c(0.05,0.95)),col="black",lty="dashed",lwd=2.5)
    abline(v=quantile(nullmodObj$Sim,c(0.025,0.975)),col="black",lty="dotted",lwd=2.5)
    mtext(as.character(date()),side=3,adj=1,line=3)
  }


