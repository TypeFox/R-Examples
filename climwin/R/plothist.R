#'Create a histogram of deltaAICc values
#'
#'Create a histogram of deltaAICc values for all fitted climate windows. Compare
#'with randomised data if provided.
#'@param dataset A dataframe containing information on all fitted climate 
#'  windows. Output from \code{\link{climatewin}}.
#'@param datasetrand A dataframe containing information on all fitted climate 
#'  windows using randomised data. Output from \code{\link{randwin}}.
#'@param histq If datasetrand is provided. The quantile of the randomised data 
#'  that will be compared to non-randomised data. Used to determine the 
#'  likelihood of finding a climate window model of a given deltaAICc value at 
#'  random.
#'@return If datasetrand is provided, plothist will return two stacked histograms
#'  to compare the deltaAICc of non-randomised and randomised data. This can 
#'  help determine the likelihood of obtaining a deltaAICc value of fitted 
#'  climate windows at random. Without datasetrand, plotall will create a single
#'  histogram of deltaAICc values for all fitted climate windows.
#'@author Liam D. Bailey and Martijn van de Pol
#'@examples
#'# Plot real and randomised data for the Mass dataset
#' 
#'data(MassOutput)
#'data(MassRand)
#' 
#'plothist(dataset = MassOutput, datasetrand = MassRand, histq = 0.95)
#'
#'# Plot deltaAICc when no randomised data is provided
#' 
#'data(MassOutput)
#' 
#'plothist(dataset = MassOutput)
#' 
#'@import ggplot2
#'@export

plothist <- function(dataset, datasetrand = NULL, histq = 0.99){
  
  if (is.null(datasetrand) == FALSE){
    
    keep2                       <- c("deltaAICc", "Randomised")
    randdata                    <- rbind(dataset[keep2], datasetrand[keep2])
    randdata$deltaAICc          <- as.numeric(randdata$deltaAICc)
    levels(randdata$Randomised) <- c("Real data", paste("Randomised data (", max(datasetrand$Repeat), "x )"))
  }

  if (is.null(datasetrand) == TRUE){
    with(dataset, {
      ggplot(dataset, aes(x = deltaAICc)) +
      geom_histogram(aes(y = 2 * ..density..), colour = "black", fill = "red", binwidth = 2, alpha = 0.5) +
      theme_classic() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(size = 0.25, colour = "black"),
            legend.position = "none",
            plot.title = element_text(size = 16))+
      ggtitle(expression(paste("Histogram of ", Delta,"AICc"))) +
      ylab("Proportion") +
      xlab(expression(paste(Delta, "AICc (compared to null model)")))
    })
  } else { 
    vline.data <- data.frame(y = as.numeric(quantile(datasetrand$deltaAICc, prob = (1 - histq))))
    with(randdata, {ggplot(randdata, aes(x = deltaAICc, fill = Randomised))+
      geom_histogram(aes(y = 2 * ..density..), colour = "black", binwidth = 2, alpha = 0.5)+
      theme_classic()+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(size = 0.25, colour = "black"),
            legend.position = "none",
            plot.title = element_text(size = 16))+
      facet_wrap(~Randomised, nrow = 2)+
      ggtitle(expression(paste("Histogram of ", Delta,"AICc")))+
      geom_vline(data = vline.data, aes(xintercept = y), linetype = "dashed", colour = "red")+
      ylab("Proportion")+
      xlab(expression(paste(Delta,"AICc (compared to null model)")))})
  }
}