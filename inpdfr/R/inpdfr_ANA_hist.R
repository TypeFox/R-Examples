#' Plot an histogram with the number of words excluding stop words
#'
#' Plot a histogram with the number of words excluding stop words using \code{\link[graphics]{hist}} function.
#' @param wordF The data.frame containing word occurrences.
#' @param mwidth The width of the plot in pixels.
#' @param mheight The height of the plot in pixels.
#' @param formatType The format for the output file ("eps", "pdf", "png", "svg", "tiff", "jpeg", "bmp").
#' @param ... Additional arguments from \code{hist} function.
#' @return NULL
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
#' getSummaryStatsHISTO(wordF = wordOccuDF)
#' @export
getSummaryStatsHISTO <- function(wordF, mwidth = 800, mheight = 800, formatType = "png", ...){

  subDir <- "RESULTS"
  dir.create(file.path(getwd(), subDir), showWarnings = FALSE)

  mFreq <- matrix(as.matrix(wordF[,2:length(wordF[1,])]), ncol = (length(wordF[1,])-1), 
    dimnames = list(as.vector(wordF[,1]), names(wordF)[2:length(wordF[1,])]))
  
  R.devices::devEval(type = formatType, name = "HISTO_numWords",
   aspectRatio = mheight / mwidth,
   scale = do.call(function(){if((mheight / mwidth) <= 1) {
     x <- max(mheight / 480, mwidth / 480)} else {
       x <- min(mheight / 480, mwidth / 480)}
     return(x)}, list())
   , path = file.path(getwd(), subDir), {
    try(graphics::hist(apply(mFreq, MARGIN = 2, FUN = sum), main = "", 
      xlab = "Number of words excluding stop words", col = grDevices::grey(0.5),...), 
      silent = TRUE)
   }
  )

  return(NULL)
}
