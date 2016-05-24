#'Plot model beta estimates
#' 
#'Create colour plots of model beta estimates. Will include quadratic and cubic
#'beta estimates where appropriate.
#'@param dataset A dataframe containing information on all fitted climate 
#' windows. Output from \code{\link{climatewin}}.
#'@param plotall Used in conjunction with function \code{\link{plotall}}. 
#' Should not be changed manually.
#'@param plotallenv Used in conjunction with function \code{\link{plotall}}.
#' Should not be changed manually.
#'@return Returns colour plots of model beta estimates. Where applicable, 2nd 
#' order coefficients (quadratic) and 3rd order coefficients (cubic) will be 
#' plotted seperately.
#'@author Liam D. Bailey and Martijn van de Pol
#'@examples
#'# Plot model beta estimates for linear models in the Mass dataset
#' 
#'data(MassOutput)
#'
#'plotbetas(dataset = MassOutput)
#' 
#'@import ggplot2
#'@import gridExtra
#'@export

plotbetas <- function(dataset, plotallenv, plotall = FALSE){
  
  with(dataset, {
    if(dataset$Function[1] == "lin" || dataset$Function[1] == "log" || dataset$Function[1] == "inv"){
      beta <-ggplot(dataset, aes(x = WindowClose, y = WindowOpen, z = ModelBeta)) +
        geom_tile(aes(fill = ModelBeta)) +
        scale_fill_gradientn(colours = c("red", "yellow", "blue"), name = "") +
        theme_classic() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line = element_line(size = 0.25, colour = "black"),
              plot.title = element_text(size = 16),
              legend.position = c(0.75,0.3)) +
        ggtitle("Beta linear") +
        ylab("Window open") +
        xlab("Window close")
      if(plotall == TRUE){
        plotallenv$beta <-beta
      } else {
        beta
      }
    } else if(dataset$Function[1] == "quad"){
      beta <-ggplot(dataset, aes(x = WindowClose, y = WindowOpen, z = ModelBeta)) +
        geom_tile(aes(fill = ModelBeta)) +
        scale_fill_gradientn(colours = c("red", "yellow", "blue"), name = "") +
        theme_classic() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line = element_line(size = 0.25, colour = "black"),
              plot.title = element_text(size = 16),
              legend.position = c(0.75,0.3)) +
        ggtitle("Beta linear") +
        ylab("Window open") +
        xlab("Window close")
      
      beta2 <- ggplot(dataset, aes(x = WindowClose, y = WindowOpen, z = ModelBetaQ)) +
        geom_tile(aes(fill = ModelBetaQ)) +
        scale_fill_gradientn(colours = c("red", "yellow", "blue"), name = "") +
        theme_classic() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line = element_line(size = 0.25, colour = "black"),
              plot.title = element_text(size = 16),
              legend.position = c(0.75, 0.3)) +
        ggtitle("Beta quadratic") +
        ylab("Window open") +
        xlab("Window close")
      if(plotall == TRUE){
        plotallenv$beta  <- beta
        plotallenv$beta2 <- beta2
      } else {
        grid.arrange(beta, beta2, nrow = 1)
      }
    } else if(dataset$Function[1] == "cub"){
      beta <-ggplot(dataset, aes(x = WindowClose, y = WindowOpen, z = ModelBeta)) +
        geom_tile(aes(fill = ModelBeta)) +
        scale_fill_gradientn(colours = c("red", "yellow", "blue"), name = "") +
        theme_classic() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line = element_line(size = 0.25, colour = "black"),
              plot.title = element_text(size = 16),
              legend.position = c(0.75,0.3)) +
        ggtitle("Beta linear") +
        ylab("Window open") +
        xlab("Window close")
      
        beta2 <- ggplot(dataset, aes(x = WindowClose, y = WindowOpen, z = ModelBetaQ)) +
          geom_tile(aes(fill = ModelBetaQ)) +
          scale_fill_gradientn(colours = c("red", "yellow", "blue"), name = "") +
          theme_classic() +
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_line(size = 0.25, colour = "black"),
                plot.title = element_text(size = 16),
                legend.position = c(0.75, 0.3)) +
          ggtitle("Beta quadratic") +
          ylab("Window open") +
          xlab("Window close")
        
        beta3 <- ggplot(dataset, aes(x = WindowClose, y = WindowOpen, z = ModelBetaC)) +
          geom_tile(aes(fill = ModelBetaC)) +
          scale_fill_gradientn(colours = c("red", "yellow", "blue"), name = "")+
          theme_classic() +
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_line(size = 0.25, colour = "black"),
                plot.title = element_text(size = 16),
                legend.position = c(0.75, 0.3)) +
          ggtitle("Beta cubic") +
          ylab("Window open") +
          xlab("Window close")
      if(plotall == TRUE){
        plotallenv$beta  <- beta
        plotallenv$beta2 <- beta2
        plotallenv$beta3 <- beta3
      } else {
        grid.arrange(beta, beta2, beta3, nrow = 1)
      }      
    } else if(dataset$Function[1] == "centre"){
      wgmean <- ggplot(dataset, aes(x = WindowClose, y = WindowOpen, z = WithinGrpMean)) +
        geom_tile(aes(fill = WithinGrpMean)) +
        scale_fill_gradientn(colours = c("red", "yellow", "blue"), name = "") +
        theme_classic() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line = element_line(size = 0.25, colour = "black"),
              plot.title = element_text(size = 16),
              legend.position = c(0.75,0.3)) +
        ggtitle("Within group mean coefficient") +
        ylab("Window open") +
        xlab("Window close")
      
      wgdev <- ggplot(dataset, aes(x = WindowClose, y = WindowOpen, z = WithinGrpDev)) +
        geom_tile(aes(fill = WithinGrpDev)) +
        scale_fill_gradientn(colours = c("red", "yellow", "blue"), name = "") +
        theme_classic() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line = element_line(size = 0.25, colour = "black"),
              plot.title = element_text(size = 16),
              legend.position = c(0.75,0.3)) +
        ggtitle("Within group deviation coefficient") +
        ylab("Window open") +
        xlab("Window close")
      
      if(plotall == TRUE){
        plotallenv$wgmean <- wgmean
        plotallenv$wgdev  <- wgdev
      } else {
        grid.arrange(wgmean, wgdev, nrow = 1)
      }
    }
  }
  )
}