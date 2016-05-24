#' Word cloud based on the word-occurrence data.frame.
#'
#' Plot a word cloud from the word-occurrence data.frame using \code{\link[wordcloud]{wordcloud}} function.
#'
#' @param wordF The data.frame containing word occurrences.
#' @param wcFormat Output format for the word cloud (deprecated, only "png").
#' @param wcminFreq Minimum word frequency for words to be ploted (see \code{\link[wordcloud]{wordcloud}}).
#' @param wcmaxWords Maximum number of words to be ploted (see \code{\link[wordcloud]{wordcloud}}).
#' @param wcRandOrder Plot words in random order (see \code{\link[wordcloud]{wordcloud}}).
#' @param wcCol Color words (see \code{\link[wordcloud]{wordcloud}}).
#' @param getPlot A vector with two logical values. If \code{plots[1]==TRUE}, a word cloud is made for each document.
#'    If \code{plots[2]==TRUE}, a word cloud is made for the combinaison of all documents.
#' @param mwidth The width of the plot in pixels.
#' @param mheight The height of the plot in pixels.
#' @param formatType The format for the output file ("eps", "pdf", "png", "svg", "tiff", "jpeg", "bmp").
#' @return NULL
#' @examples
#' \dontrun{
#' makeWordcloud(wordF = myDF)
#' }
#' @export
makeWordcloud <- function(wordF, wcFormat = "png", wcminFreq = 3, wcmaxWords = Inf,
  wcRandOrder = FALSE, wcCol = RColorBrewer::brewer.pal(8, "Dark2"),
  getPlot = c(TRUE, TRUE), mwidth = 1000, mheight = 1000, formatType = "png"){

  subDir <- "RESULTS"
  dir.create(file.path(getwd(), subDir), showWarnings = FALSE)

  lapply(2:length(wordF), function(i){
    if (wcFormat == "png" && getPlot[1] == TRUE){
      R.devices::devEval(type = formatType, name = paste0("WORDCLOUD_",names(wordF)[i]),
       aspectRatio = mheight / mwidth,
       scale = do.call(function(){if((mheight / mwidth) <= 1) {
         x <- max(mheight / 480, mwidth / 480)} else {
           x <- min(mheight / 480, mwidth / 480)}
         return(x)}, list())
       , path = file.path(getwd(), subDir), {
        try(wordcloud::wordcloud(words = wordF[,1], freq = wordF[,i], 
          min.freq = min(wcminFreq, max(wordF[,i])), max.words = wcmaxWords,
          random.order = wcRandOrder, col = wcCol), silent = TRUE)
       }
      )
    }
  })
  if(length(wordF[1,]) > 2){
    wordFreqMerged <- rowSums(wordF[,2:length(wordF[1,])])
  }else{
    wordFreqMerged <- wordF[,2]
  }
  wordFreqMerged <- data.frame(wordF[,1], wordFreqMerged)
  if (getPlot[2] == TRUE){
    R.devices::devEval(type = formatType, name = "WORDCLOUD_ALL",
     aspectRatio = mheight / mwidth,
     scale = do.call(function(){if((mheight / mwidth) <= 1) {
       x <- max(mheight / 480, mwidth / 480)} else {
         x <- min(mheight / 480, mwidth / 480)}
       return(x)}, list())
     , path = file.path(getwd(), subDir), {
      try(wordcloud::wordcloud(words = wordFreqMerged[,1], freq = wordFreqMerged[,2], 
        min.freq = min(wcminFreq, max(wordFreqMerged[,2])), max.words = wcmaxWords,
        random.order = wcRandOrder, col = wcCol), silent = TRUE)
     }
    )
  }

  return(NULL)
}
