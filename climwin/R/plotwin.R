#'Plot the opening and closing point of best climate windows
#'
#'Visualise the opening and closing point for a subset of best climate windows.
#'@param dataset A dataframe containing information on all fitted climate
#'  windows. Output from \code{\link{climatewin}}.
#'@param cw Cumulative model weight used to subset the group of best models.
#'@return Creates two boxplots showing the opening and closing point for a subset
#'  of best climate windows. Best climate windows make up the
#'  cumulative model weight equivalent to the value of cw.
#'@author Liam D. Bailey and Martijn van de Pol
#'@examples
#'# View window limits for climate windows in the top 95% of model weights.
#' 
#'data(MassOutput)
#' 
#'plotwin(dataset = MassOutput, cw = 0.95)
#' 
#'@import ggplot2
#'@import reshape
#'@export

plotwin <- function(dataset, cw = 0.95){
  
  #Order models by weight#
  dataset    <- dataset[order(-dataset$ModWeight), ]
  dataset$cw <- as.numeric(cumsum(dataset$ModWeight) <= cw)
  datasetcw  <- subset(dataset, cw == 1)
  
  keep=c("closest", "WindowClose", "WindowOpen")
  
  datasetcw                  <- datasetcw[keep]
  datasetcw                  <- melt(datasetcw, id = "closest")
  datasetcw$variable         <- factor(datasetcw$variable, levels = c("WindowOpen", "WindowClose"))
  levels(datasetcw$variable) <- c("Window Open", "Window Close")
  
  p_meds <- data.frame(variable = levels(datasetcw$variable), value = as.numeric(tapply(datasetcw$value, datasetcw$variable, median)))
  
  with(datasetcw, {
    ggplot(datasetcw, aes(x = variable, y = value))+
      scale_y_continuous(limits = c(min(dataset$WindowClose), max(dataset$WindowClose)))+
      geom_boxplot(width = 0.5)+
      geom_text(data = p_meds, aes(x = variable, y = value, label = value),
                size = 5, vjust = -1.9) +
      coord_flip()+
      theme_classic()+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(size = 0.25, colour = "black"),
            axis.text.y = element_text(angle = 90, hjust = 0.5,size = 10),
            plot.title = element_text(size = 16))+
      ggtitle(paste("Climate window range for top \n", (cw*100), "% of model weights"))+
      xlab("")+
      ylab("Climate window")
  })
}