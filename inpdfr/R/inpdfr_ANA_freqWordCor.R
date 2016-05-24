#' Returns most frequent words.
#'
#' Returns most frequent words and plots their frequencies per document.
#'
#' @param wordF The data.frame containing word occurrences.
#' @param numWords The number of words to be returned.
#' @param getPlot If \code{TRUE}, save a scatter plot in the RESULTS directory.
#' @param mwidth The width of the plot in pixels.
#' @param mheight The height of the plot in pixels.
#' @param formatType The format for the output file ("eps", "pdf", "png", "svg", "tiff", "jpeg", "bmp").
#' @return The \code{numWords} most frequent words.
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
#' getMostFreqWord(wordF = wordOccuDF, numWords = 5)
#' @export
getMostFreqWord <- function(wordF, numWords, getPlot = TRUE, mwidth = 1024, mheight = 800, 
  formatType = "png"){

  subDir <- "RESULTS"
  dir.create(file.path(getwd(), subDir), showWarnings = FALSE)

  mostFreqWords <- NULL
  if(is.numeric(numWords)){
    if(numWords > nrow(wordF)){numWords <- nrow(wordF)}
    mostFreqWords <- wordF[1:numWords, 1]
    mostFreqWords <- as.character(mostFreqWords)
  }
  mXwords <- wordF[1:numWords, 2:ncol(wordF)]

  if(getPlot == TRUE){
    R.devices::devEval(type = formatType, name = paste0('MostFreqWords_',numWords),
     aspectRatio = mheight / mwidth,
     scale = do.call(function(){if((mheight / mwidth) <= 1) {
       x <- max(mheight / 480, mwidth / 480)} else {
         x <- min(mheight / 480, mwidth / 480)}
       return(x)}, list())
     , path = file.path(getwd(), subDir), {
      graphics::plot(1, type = "n", xlim = c(1, ncol(mXwords)), ylim = c(0, max(mXwords)), 
        ylab = "Occurrences", xlab = "", axes = FALSE)
      graphics::axis(1, at = 1:ncol(mXwords), labels = names(wordF[2:ncol(wordF)]), las = 2)
      graphics::axis(2)
      mycol <- 1:numWords
      sapply(1:nrow(mXwords), function(i){graphics::points(x = 1:ncol(mXwords), 
        y = as.vector(mXwords[i,]), type = "o", col = mycol[i], lwd = 2)})
      graphics::legend("topright", legend = as.character(wordF[1:numWords, 1]), lty = 1, 
        col = mycol, lwd = 2)
     }
    )
  }
  return(mostFreqWords)
}

#' Test for correlation between the most frequent words.
#'
#' @param wordF The data.frame containing word occurrences.
#' @param numWords The number of words to be returned.
#' @param getPlot A vector with two logical values. If \code{plots[1]==TRUE},
#'    an image of the correlation matrix is saved in the RESULTS directory.
#'    If \code{plots[2]==TRUE}, the image of the p-value matrix associated
#'    with the correlation is saved in the RESULTS directory.
#' @param getTextSink If \code{TRUE}, save the correlation matrix and the
#'    associated p-values in a text file in the RESULTS directory.
#' @param mwidth The width of the plot in pixels.
#' @param mheight The height of the plot in pixels.
#' @param formatType The format for the output file ("eps", "pdf", "png", "svg", "tiff", "jpeg", "bmp").
#' @return A list with the correlation matrix and the p-value matrix.
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
#' getMostFreqWordCor(wordF = wordOccuDF, numWords = 5)
#' @export
getMostFreqWordCor <- function(wordF, numWords, getPlot = c(TRUE, TRUE), getTextSink = TRUE,
  mwidth = 1024, mheight = 1024, formatType = "png"){ # correlation between words

  subDir <- "RESULTS"
  dir.create(file.path(getwd(), subDir), showWarnings = FALSE)

  if(numWords > nrow(wordF)){numWords <- nrow(wordF)}
  M <- wordF[,2:ncol(wordF)]
  matCOR <- matrix(NA, ncol = numWords, nrow = numWords)
  for (i in 1:numWords){
    for (j in 1:numWords){
      matCOR[i,j] <- stats::cor(as.vector(unlist(M[i,])), as.vector(unlist(M[j,])))
    }
  }
  colnames(matCOR) <- wordF[1:numWords,1]
  rownames(matCOR) <- wordF[1:numWords,1]

  matCORtest <- matrix(NA, ncol = numWords, nrow = numWords)
  for (i in 1:numWords){
    for (j in 1:numWords){
      matCORtest[i,j] <- stats::cor.test(as.vector(unlist(M[i,])), as.vector(unlist(M[j,])))$p.value
    }
  }
  colnames(matCORtest) <- wordF[1:numWords,1]
  rownames(matCORtest) <- wordF[1:numWords,1]
  # matCORtest[matCORtest>0.05]<-NA

  if(getPlot[2] == TRUE){
    R.devices::devEval(type = formatType, name = paste0('MostFreqWordsCorPvalue_',numWords),
     aspectRatio = mheight / mwidth,
     scale = do.call(function(){if((mheight / mwidth) <= 1) {
       x <- max(mheight / 480, mwidth / 480)} else {
         x <- min(mheight / 480, mwidth / 480)}
       return(x)}, list())
     , path = file.path(getwd(), subDir), {
      graphics::par(mar = c(7, 7, 1, 1))
      graphics::image(matCORtest, axes = FALSE, col = grDevices::heat.colors(5))
      graphics::axis(1, at = seq(0, 1, length = numWords), labels = colnames(matCOR), las = 2)
      graphics::axis(2, at = seq(0, 1, length = numWords), labels = colnames(matCOR), las = 1)
     }
    )
  }

  matCorSign <- as.data.frame(matCOR)
  matCorSign[matCOR <= 1] <- "****"
  matCorSign[matCOR < 0.999] <- "***"
  matCorSign[matCOR < 0.75] <- "**"
  matCorSign[matCOR < 0.50] <- "*"
  matCorSign[matCOR < 0.25] <- "."
  matCorSign[matCOR < 0.10] <- ""
  matCorSign[matCOR < (-0.15)] <- ""
  matCorSign[matCOR < (-0.25)] <- "(*)"
  matCorSign[matCOR < (-0.50)] <- "(**)"
  matCorSign[matCOR < (-0.75)] <- "(***)"

  if(getTextSink == TRUE){
    sink(paste0('RESULTS/MostFreqWordsCor_', numWords, '.txt'))
    cat('\n#######################\n### RAW             ###\n#######################\n')
    try(print(matCOR), silent = TRUE)
    cat('\n#######################\n### SIGN            ###\n#######################\n')
    cat('# -1;-0.75		(***)\n# -0.75;-0.5	(**)\n# -0.5;-0.25	(*)\n# -0.25;-0.10	(.)\n# -0.10;0.10	\n# 0.10;0.25		.\n# 0.25;0.5		*\n# 0.5;0.75		**\n# 0.75;0.999	***\n# 0.999;1		****\n\n')
    try(print(matCorSign), silent = TRUE)
    cat('\n#######################\n### PVALUE          ###\n#######################\n')
    try(print(matCORtest), silent = TRUE)
    sink()
  }

  if(getPlot[1] == TRUE){
    R.devices::devEval(type = formatType, name = paste0('MostFreqWordsCor_',numWords),
     aspectRatio = mheight / mwidth,
     scale = do.call(function(){if((mheight / mwidth) <= 1) {
       x <- max(mheight / 480, mwidth / 480)} else {
         x <- min(mheight / 480, mwidth / 480)}
       return(x)}, list())
     , path = file.path(getwd(), subDir), {
      graphics::par(mar = c(7, 7, 1, 1))
      graphics::image(abs(matCOR), axes = FALSE, col = rev(grDevices::heat.colors(5)))
      graphics::axis(1, at = seq(0, 1, length = numWords), labels = colnames(matCOR), las = 2)
      graphics::axis(2, at = seq(0, 1, length = numWords), labels = colnames(matCOR), las = 1)
     }
    )
  }

  return(list(cor = matCOR, pval = matCORtest))
}

#' Returns most frequent words
#'
#' @param wordF The data.frame containing word occurrences.
#' @param occuWords The minimum number of occurrences for words to be returned.
#' @return A vector with most frequent words.
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
#' getXFreqWord(wordF = wordOccuDF, occuWords = 5)
#' @export
getXFreqWord <- function(wordF, occuWords){
  xFreqWords <- NULL
  if(is.numeric(occuWords)){
    datasetSum <- apply(wordF[,2:ncol(wordF)], MARGIN = 1, FUN = sum)
    xFreqWords <- wordF[,1][datasetSum >= occuWords]
    xFreqWords <- as.character(xFreqWords)
  }
  return(xFreqWords)
}




