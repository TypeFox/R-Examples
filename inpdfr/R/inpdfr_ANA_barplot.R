#' Perform a barplot with the number of unique words per document
#'
#' Perform a barplot with the number of unique words per document using \code{\link[graphics]{barplot}} function.
#' @param wordF The data.frame containing word occurrences.
#' @param getPlot If \code{TRUE}, save the bar plot in the RESULTS directory.
#' @param mwidth The width of the plot in pixels.
#' @param mheight The height of the plot in pixels.
#' @param formatType The format for the output file ("eps", "pdf", "png", "svg", "tiff", "jpeg", "bmp").
#' @param ... Additional arguments from \code{barplot} function.
#' @return The number of unique words per document.
#' @examples
#' data("loremIpsum")
#' loremIpsum01 <- loremIpsum[1:100]
#' loremIpsum02 <- loremIpsum[101:200]
#' loremIpsum03 <- loremIpsum[201:300]
#' loremIpsum04 <- loremIpsum[301:400]
#' loremIpsum05 <- loremIpsum[401:500]
#' subDir <- "RESULTS"
#' dir.create(file.path(getwd(), subDir), showWarnings = FALSE)
#' write(x = loremIpsum01, file = "RESULTS/loremIpsum01.txt")
#' write(x = loremIpsum02, file = "RESULTS/loremIpsum02.txt")
#' write(x = loremIpsum03, file = "RESULTS/loremIpsum03.txt")
#' write(x = loremIpsum04, file = "RESULTS/loremIpsum04.txt")
#' write(x = loremIpsum05, file = "RESULTS/loremIpsum05.txt")
#' wordOccuDF <- getwordOccuDF(mywd = paste0(getwd(), "/RESULTS"), 
#'   excludeSW = FALSE)
#' file.remove(list.files(pattern = "loremIpsum"))
#' getSummaryStatsBARPLOT(wordF = wordOccuDF)
#' @export
getSummaryStatsBARPLOT <- function(wordF, getPlot = TRUE, mwidth = 480, 
  mheight = 480, formatType = "png", ...){
  ## create RESULTS folder
  subDir <- "RESULTS"
  dir.create(file.path(getwd(), subDir), showWarnings = FALSE)
  ## make barplot
  getBarPlot <- function(mFreq, seqBar){
    mFreqSup1 <- mFreq
    mFreqSup1[mFreqSup1 <= max(mFreq)/seqBar[1]] <- 0
    mFreqSup2 <- mFreq
    mFreqSup2[mFreqSup2 <= max(mFreq)/seqBar[2]] <- 0
    mFreqSup3 <- mFreq
    mFreqSup3[mFreqSup3 <= max(mFreq)/seqBar[3]] <- 0
    mFreqSup3Count <- sapply(1:ncol(mFreqSup3), function(i){length(mFreqSup3[,i][mFreqSup3[,i] != 0])})
    mFreqSup2Count <- sapply(1:ncol(mFreqSup2), function(i){length(mFreqSup2[,i][mFreqSup2[,i] != 0])})
    mFreqSup1Count <- sapply(1:ncol(mFreqSup1), function(i){length(mFreqSup1[,i][mFreqSup1[,i] != 0])})
    mFreqSup0Count <- sapply(1:ncol(mFreq), function(i){length(mFreq[,i][mFreq[,i] != 0])})
    barPlot1 <- rbind(
      mFreqSup3Count,
      mFreqSup2Count - mFreqSup3Count,
      mFreqSup1Count - mFreqSup2Count,
      mFreqSup0Count - mFreqSup1Count,
      ...
    )
    return(barPlot1)
  }
  mFreq <- matrix(as.matrix(wordF[,2:length(wordF[1,])]),
    ncol = (length(wordF[1,]) - 1),
    dimnames = list(as.vector(wordF[,1]),
      names(wordF)[2:length(wordF[1,])]
    )
  )
  barPlot1 <- getBarPlot(mFreq = mFreq, seqBar = c(8, 3, 2))
  ## save barplot
  if(getPlot == TRUE){
    R.devices::devEval(type = formatType, name = "BARPLOT_numWords",
      aspectRatio = mheight / mwidth,
      scale = do.call(function(){if((mheight / mwidth) <= 1) {
        x <- max(mheight / 480, mwidth / 480)} else {
        x <- min(mheight / 480, mwidth / 480)}
        return(x)}, list())
        , path = file.path(getwd(), subDir), {
        ## graph
        graphics::par(mar = c(10, 4, 1, 1))
        try(graphics::barplot(barPlot1, las = 2, ylab = "Number of unique words",
          col = c(grDevices::grey(0.3), grDevices::grey(0.5), grDevices::grey(0.7), grDevices::grey(0.9)),
          names = colnames(mFreq)), silent = TRUE
        )
        try(graphics::legend("topright", bg = "white",
          fill = c(grDevices::grey(0.3), grDevices::grey(0.5), grDevices::grey(0.7), grDevices::grey(0.9)),
          legend = c(paste0("Words used more than: ", floor(max(mFreq) / 2)),
            paste0("Words used more than: ", floor(max(mFreq) / 3)),
            paste0("Words used more than: ", floor(max(mFreq) / 8)),
            paste0("Words used less than or equal to: ",
              floor(max(mFreq) / 8)
            )
          )), silent = TRUE
        )
        ## end graph
    })
  }
  ## return raw data
  return(apply(barPlot1, MARGIN = 2, FUN = sum, na.rm = TRUE))
}
