#' Truncate the word-occurrence data.frame.
#'
#' @param wordF The data.frame containing word occurrences.
#' @param maxWords The maximum number of words in the data.frame.
#' @return The data.frame containing word occurrences.
#' @examples
#' \dontrun{
#' truncNumWords(wordF = myWordOccurrenceDF, maxWords = 50)
#' }
#' @export
truncNumWords <- function(wordF, maxWords){
  try(if((!is.na(maxWords) || maxWords != "")==TRUE){
    for(i in 1:length(wordF)){
      numWords <- length(wordF[[i]][[1]][,1])
      if(numWords<maxWords){maxWords_i <- numWords}else{maxWords_i <- maxWords}
      wordF[[i]][[1]] <- wordF[[i]][[1]][1:maxWords_i,]
    }
  }, silent = TRUE)
  return(wordF)
}

#' Merge word-occurrence data.frames into a single data.frame.
#'
#' @param wordF The data.frame containing word occurrences.
#' @return A single word-occurrrence data.frame with each column corresponding to a text file.
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
#' wordOccuFreq <- getTXT(myTXTs = list.files(path = paste0(getwd(), 
#'   "/RESULTS/"), pattern = "loremIpsum", full.names = TRUE))
#' wordOccuDF <- mergeWordFreq(wordF = wordOccuFreq)
#' file.remove(list.files(pattern = "loremIpsum"))
#' @export
mergeWordFreq <- function(wordF){
  fileNames <- sapply(wordF, "[[", 2)
  fileWordFreq <- lapply(wordF, "[[", 1)
  words <- unique(as.vector(unique(unlist(sapply(fileWordFreq, function(i){unique(i$word)})))))
  mydbWords <- data.frame(word = words)
  for(k in 1:length(fileWordFreq)){
    mydbWords <- try(merge(mydbWords, fileWordFreq[[k]][,c(1, 3)], by.x = 1, by.y = 2, all = TRUE, suffixes = c(paste0(".x", k), paste0(".y", k))), silent = TRUE)
  }
  colnames(mydbWords) <- c("word", fileNames)
  mydbWords[is.na(mydbWords)] <- 0

  if(length(fileNames)>1){
    ### stem over the incidence-matrix
    dd <- NULL
    # Stem words
    dd$stem <- SnowballC::wordStem(as.character(mydbWords[,1]), language = "english") ###
    dd$word <- as.character(mydbWords[,1])
    dd$freq <- mydbWords[,2:ncol(mydbWords)]

    agg_freq <- sapply(1:ncol(dd$freq), function(i){stats::aggregate(freq[,i] ~ stem, data = dd, sum)}) ### a optimiser
    agg_freq <- cbind(data.frame(unlist(agg_freq[1,1])), data.frame(agg_freq[2,]))
    agg_word <- stats::aggregate(word ~ stem, data = dd, function(x) x[1]) ###
    dd <- cbind(word = agg_word[,2], agg_freq[,2:ncol(agg_freq)]) ###
    names(dd) <- names(mydbWords)
    mydbWords <- dd
    ### end stem
  }
  if(length(fileNames)>1){
    mydbWords <- mydbWords[order(apply(mydbWords[,2:ncol(mydbWords)], MARGIN = 1, FUN = sum), decreasing = T), ] # order matrix by incidence
  }
  return(mydbWords)
}
