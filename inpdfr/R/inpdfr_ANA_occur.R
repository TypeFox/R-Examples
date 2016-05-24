#' Plot a scatter plot with the proportion of documents using similar words.
#'
#' @param wordF The data.frame containing word occurrences.
#' @param getPlot If \code{TRUE}, save the scatter plot in the RESULTS directory.
#' @param mwidth The width of the plot in pixels.
#' @param mheight The height of the plot in pixels.
#' @param formatType The format for the output file ("eps", "pdf", "png", "svg", "tiff", "jpeg", "bmp").
#' @return A data.frame containing the proportion of documents and the number of similar words.
#' @examples
#' \dontrun{
#' getSummaryStatsOCCUR(wordF = myDF)
#' }
#' @export
getSummaryStatsOCCUR <- function(wordF, getPlot = TRUE, mwidth = 800, mheight = 800, 
  formatType = "png"){

  subDir <- "RESULTS"
  dir.create(file.path(getwd(), subDir), showWarnings = FALSE)

  ncores <- parallel::detectCores()
  if(nrow(wordF) < ncores){ncores <- nrow(wordF)}
  cl<-parallel::makeCluster(ncores)
  # parallel::clusterExport(cl=cl, varlist=c("wordF")) ### for testing purposes
  on.exit(parallel::stopCluster(cl))
  tableP1 <- parallel::parSapply(cl, 1:nrow(wordF), function(i){
    length(wordF[i, 2:ncol(wordF)][wordF[i, 2:ncol(wordF)] != 0]) / (ncol(wordF) - 1)
  })
  tableP2 <- table(tableP1)
  dfTableP <- data.frame(numW = tableP2, breaks = cut(as.numeric(names(tableP2)), 
    breaks = seq(0, 1, 0.1)))
  aggregDfTableP <- stats::aggregate(dfTableP[,2] ~ dfTableP[,3], FUN = sum)

  if(getPlot == TRUE){
    R.devices::devEval(type = formatType, name = "PLOT_occurrence",
     aspectRatio = mheight / mwidth,
     scale = do.call(function(){if((mheight / mwidth) <= 1) {
       x <- max(mheight / 480, mwidth / 480)} else {
         x <- min(mheight / 480, mwidth / 480)}
       return(x)}, list())
     , path = file.path(getwd(), subDir), {
      graphics::par(mar = c(5, 4, 2, 2))
      try(graphics::plot(aggregDfTableP, las = 2, log = "y", ylab = "Number of words", xlab = ""), 
        silent = TRUE)
     }
    )
  }

  return(aggregDfTableP)
}
