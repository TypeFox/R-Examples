#'Visualise climate window data
#'
#'Creates a panel of plots to help visualise climate window data.
#'@param dataset A dataframe containing information on all fitted climate 
#'  windows. Output from \code{\link{climatewin}}.
#'@param datasetrand A dataframe containing information on all fitted climate 
#'  windows using randomised data. Output from \code{\link{randwin}}.
#'@param bestmodel A model object. The strongest climate window model. Returned 
#'  from \code{\link{singlewin}} or \code{\link{climatewin}}.
#'@param bestmodeldata A dataframe containing the biological and climate data
#'  used to fit the strongest climate window model. Output from
#'  \code{\link{singlewin}} or \code{\link{climatewin}}.
#'@param cw1,cw2,cw3 Cumulative weight levels used to visualise model weight 
#'  distribution. See \code{\link{plotweights}} for more detail.
#'@param histq If datasetrand is provided. The quantile of the randomised data 
#'  to be compared with non-randomised data. Can be used to determine the 
#'  likelihood of finding a climate window model of a given deltaAICc value by
#'  chance.
#'@param title Title of the plot panel.
#'@return Will return a panel of 6-8 plots:
#'  
#'  \itemize{
#'  \item DeltaAICc: A colour plot of model deltaAICc values (larger
#'  negative values indicate stronger models). DeltaAICc is the difference
#'  between AICc of each climate window model and the baseline model containing
#'  no climate data.
#'  
#'  \item Model weight: A plot showing the distribution of cumulative
#'  model weights. Gradient levels determined by parameters cw1, cw2 and cw3.
#'  Darker areas have a higher chance of containing the best climate window.
#'  
#'  \item Model betas: A colour plot of model beta estimates. Where applicable,
#'  2nd order coefficients (quadratic) and 3rd order coefficients (cubic) will
#'  be plotted separately.
#'  
#'  \item Histogram(s): If datasetrand is provided, plotall will create two 
#'  stacked histograms to compare the deltaAICc of non-randomised and randomised
#'  data. This can help determine the likelihood of obtaining a deltaAICc value 
#'  for a fitted climate window model at random. Without datasetrand, plotall
#'  will create a single histogram of deltaAICc values for all fitted climate 
#'  windows.
#'  
#'  \item Boxplots: Two boxplots showing the opening and closing day for a 
#'  subset of best climate windows. Best climate windows make up the
#'  cumulative model weight equivalent to the largest value of cw1, cw2 and cw3.
#'  Values above boxplots represent the median values.
#'  
#'  \item Best Model: If bestmodel and bestmodeldata are provided, plotall will 
#'  create a scatterplot to show the fit of the best model through the data. }
#'  
#'@author Liam D. Bailey and Martijn van de Pol
#'@examples
#'
#'# Visualise a fixed climate window generated for dataframes Mass and MassClimate
#'
#'data(MassOutput)
#'data(Mass)
#'data(MassClimate)
#'
#'single <- singlewin(xvar = list(Temp = MassClimate$Temp), 
#'                    cdate = MassClimate$Date, bdate = Mass$Date, 
#'                    baseline = lm(Mass ~ 1, data = Mass), 
#'                    furthest = 72, closest = 15, 
#'                    stat = "mean", func = "lin", 
#'                    type = "fixed", cutoff.day = 20, cutoff.month = 5, 
#'                    cmissing = FALSE, cinterval = "day")
#'            
#'plotall(dataset = MassOutput, bestmodel = single$bestmodel, 
#'        bestmodeldata = single$bestmodeldata,
#'        cw1 = 0.95, cw2 = 0.5, cw3 = 0.25, histq = 0.99, title = "Mass")
#'         
#'          
#'@import gridExtra
#'@import ggplot2
#'@export

plotall <- function(dataset, datasetrand = NULL,
                    bestmodel = NULL, bestmodeldata = NULL,
                    cw1 = 0.95, cw2 = 0.5, cw3 = 0.25, histq = 0.99,
                    title = NULL){
  
  a       <- c(cw1, cw2, cw3)
  b       <- a[order (-a)]
  cwa     <- b[1]
  cwb     <- b[2]
  cwc     <- b[3]
  plotenv <- environment()
  
  plotbetas(dataset = dataset, plotall = TRUE, plotallenv = plotenv)
  
  delta  <- plotdelta(dataset = dataset)
  
  cw     <- plotweights(dataset = dataset, cw1 = cwa, cw2 = cwb, cw3 = cwc)
  
  window <- plotwin(dataset = dataset, cw = cwa)
  
  hist   <- plothist(dataset = dataset, datasetrand = datasetrand, histq = histq)
  if(is.null(bestmodel) == FALSE && is.null(bestmodeldata) == FALSE){
  best  <- plotbest(dataset = dataset, bestmodel = bestmodel, bestmodeldata = bestmodeldata)
  
  if (dataset$Function[1] == "lin"){
    gridExtra::grid.arrange(delta, cw, plotenv$beta, hist, window, best, nrow = 2, ncol = 3, top = paste(title))
  } else if (dataset$Function[1] == "quad"){
    gridExtra::grid.arrange(delta, cw, plotenv$beta, plotenv$beta2, hist, window, best, nrow = 2, ncol = 4, top = paste(title))
  } else if(dataset$Function[1] == "cub"){
    gridExtra::grid.arrange(delta, cw, plotenv$beta, plotenv$beta2, hist, window, best, plotenv$beta3, nrow = 2, ncol = 4, top = paste(title))
  } else if(dataset$Function[1] == "centre"){
    gridExtra::grid.arrange(delta, cw, plotenv$wgmean, plotenv$wgdev, hist, window, best, nrow = 2, ncol = 4, top = paste(title))
  } else {
    gridExtra::grid.arrange(plotenv$beta, delta, cw, hist, window, best, nrow = 2, top = paste(title))
  }
  } else {
    if (dataset$Function[1] == "lin"){
      gridExtra::grid.arrange(delta, cw, plotenv$beta, hist, window, nrow = 2, ncol = 3, top = paste(title))
    } else if (dataset$Function[1] == "quad"){
      gridExtra::grid.arrange(delta, cw, plotenv$beta, plotenv$beta2, hist, window, nrow = 2, ncol = 4, top = paste(title))
    } else if(dataset$Function[1] == "cub"){
      gridExtra::grid.arrange(delta, cw, plotenv$beta, plotenv$beta2, hist, window, plotenv$beta3, nrow = 2, ncol = 4, top = paste(title))
    } else if(dataset$Function[1] == "centre"){
      gridExtra::grid.arrange(delta, cw, plotenv$wgmean, plotenv$wgdev, hist, window, nrow = 2, ncol = 4, top = paste(title))
    } else {
      gridExtra::grid.arrange(plotenv$beta, delta, cw, hist, window, nrow = 2, top = paste(title))
    } 
  }
}