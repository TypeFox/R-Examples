#' A quick way to obtain the word-occurrence data.frame from a set of documents.
#'
#' @param mywd A character variable containing the working directory.
#' @param language The language used ("French", "English", "Spanish").
#' @param excludeSW A logical to exclude stop words.
#' @return A single word-occurrrence data.frame.
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
#' @export
getwordOccuDF <- function(mywd, language = "English", excludeSW = TRUE){
  listFilesExt <- getListFiles(mywd)
  if (length(listFilesExt$pdf) > 0){
    wordFreqPDF <- getPDF(myPDFs = listFilesExt$pdf)
  } else {wordFreqPDF <- NULL}
  wordFreqTXT <- getTXT(myTXTs = listFilesExt$txt)
  wordFreq <- append(wordFreqPDF, wordFreqTXT)
  if (excludeSW == TRUE){
    wordFreq <- excludeStopWords(wordF = wordFreq, lang = language)
  }
  wordFreq <- truncNumWords(maxWords = Inf, wordF = wordFreq)
  mergedD <- mergeWordFreq(wordF = wordFreq)
  return(mergedD)
}

#' A quick way to compute a set of analysis from the word-occurrence data.frame.
#'
#' @param dataset A single word-occurrrence data.frame.
#' @param wcloud A logical to for word cloud analysis.
#' @param sumStats A logical to for summary statistics analysis.
#' @param freqW A logical to for word frequency analysis.
#' @param corA A logical to for correspondence analysis.
#' @param clust A logical to for cluster analysis.
#' @param metacom A logical to for metacommunity analysis.
#' @return A set of analyses available from the \code{inpdfr} package.
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
#' getAllAnalysis(dataset = wordOccuDF, wcloud = FALSE, sumStats = FALSE)
#' @export
getAllAnalysis <- function(dataset, wcloud = TRUE, sumStats = TRUE, freqW = TRUE,
  corA = TRUE, clust = TRUE, metacom = TRUE){
  if(wcloud == TRUE){
    makeWordcloud(wordF = dataset, wcminFreq = 50, wcmaxWords = Inf, 
      wcRandOrder = FALSE, wcCol = RColorBrewer::brewer.pal(8, "Dark2"), 
      getPlot = c(FALSE, TRUE))
  }
  if(sumStats == TRUE){
    getSummaryStatsBARPLOT(wordF = dataset)
    getSummaryStatsHISTO(wordF = dataset)
    getSummaryStatsOCCUR(wordF = dataset)
  }
  if(freqW == TRUE){
    getMostFreqWord(wordF = dataset, numWords = 5)
    getMostFreqWord(wordF = dataset, numWords = 50)
    getMostFreqWord(wordF = dataset, numWords = 100)
    getXFreqWord(wordF = dataset, 50)
  }
  if(corA == TRUE){
    doCA(wordF = dataset)
  }
  if(clust == TRUE){
    doCluster(wordF = dataset, myMethod = "ward.D2", gp = FALSE, nbGp = 5)
    doKmeansClust(wordF = dataset, nbClust = 4, nbIter = 10, algo = "Hartigan-Wong")
  }
  if(metacom == TRUE){
    doMetacomEntropart(wordF = dataset)
    doMetacomMetacom(wordF = dataset, numSim = 10, limit = "Inf")
  }
}
