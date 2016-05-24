#' Calculate percent variance of eigenvalues for plot_shiny.mfpca()
#'
#' Internal method that calculates percent variance of eigenvalues for specified level (1, 2, or total) for 
#' plot_shiny.mfpca(). The desired level is passed in as an argument (level = 12 for total) and a list of percent
#' variances is returned.
#' 
#' @param level numeric, 1 or 2 for levels 1 or 2, respectively, 12 to calculate total variance.
#' @param plotObj the mfpca object plotted in the plot_shiny.mfpca() function.
#' @return a list of numbers that indicate percent variance for selected level.
#' @author Julia Wrobel \email{jw3134@@cumc.columbia.edu}
#' 
varPercent <- function(level, plotObj){
  if (level == 12){
    eigenvalues = c(plotObj$evalues$level1, plotObj$evalues$level2)
  }else{eigenvalues = plotObj$evalues[[level]]}
  
  percent <- lapply(eigenvalues, function(i){100*round(i/sum(eigenvalues), 3)})
  return(percent)
}



