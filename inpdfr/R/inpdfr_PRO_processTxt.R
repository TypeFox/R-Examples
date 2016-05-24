#' Extract text from txt files and pre-process content.
#'
#' @param filetxt A character containing the name of a txt file.
#' @param encodingIn Encoding of the text file (default = "UTF-8").
#' @param encodingOut Encoding of the text extracted (default = "UTF-8").
#' @return A character vector with the content of the pre-process txt file (one element per line).
#' @examples
#' data("loremIpsum")
#' subDir <- "RESULTS"
#' dir.create(file.path(getwd(), subDir), showWarnings = FALSE)
#' write(x = loremIpsum, file = "RESULTS/loremIpsum.txt")
#' preProcTxt(filetxt = paste0(getwd(), "/RESULTS/loremIpsum.txt"))
#' file.remove(list.files(pattern = "loremIpsum"))
#' @export
preProcTxt <- function(filetxt, encodingIn = "UTF-8", encodingOut = "UTF-8"){
  zz <- file(filetxt, 'r', encoding = encodingIn)
  txt <- readChar(zz, file.info(filetxt)$size)
  close(zz)
  Encoding(txt) <- encodingOut
  txt <- strsplit(txt, split = "\\.")[[1]]
  txt <- tm::stripWhitespace(txt) # remove extra-spaces
  txt <- tolower(txt) # lowercase txt
  txt <- gsub("\\d", "", txt) # remove numbers
  txt <- gsub("'", " ", txt) # remove "'" and replace by space
  txt <- gsub("-", " ", txt) # remove "-" and replace by space
  txt <- gsub("\f", "", txt) # remove "\f"
  txt <- gsub("([a-z])\\1\\1\\1", "", txt) # remove three consecutive same letters
  txt <- gsub("([a-z])\\1\\1\\1", "", txt) # again, if previous step generated new patterns
  return(txt)
}

#' Prossess vectors containing words into a data.frame of word occurrences.
#'
#' @param txt A vector containing text.
#' @param minword An integer specifying the minimum number of letters per word into the returned data.frame.
#' @param maxword An integer to specifying the maximum number of letters per word into the returned data.frame.
#' @param minFreqWord An integer specifying the minimum word frequency into the returned data.frame.
#' @return A data.frame (freq = occurrences, stem = stem words, word = words), sorted by word occurrences.
#' @examples
#' \dontrun{
#' postProcTxt(txt = preProcTxt(filetxt = "loremIpsum.txt"))
#' }
#' data("loremIpsum")
#' subDir <- "RESULTS"
#' dir.create(file.path(getwd(), subDir), showWarnings = FALSE)
#' write(x = loremIpsum, file = "RESULTS/loremIpsum.txt")
#' preProcTxt(filetxt = paste0(getwd(), "/RESULTS/loremIpsum.txt"))
#' postProcTxt(txt = preProcTxt(filetxt = paste0(getwd(), "/RESULTS/loremIpsum.txt")))
#' file.remove(list.files(pattern = "loremIpsum"))
#' @export
postProcTxt <- function(txt, minword = 1, maxword = 20, minFreqWord = 1){
  txtMerged <- paste(txt, collapse = ' ')
  txtMerged <- gsub("[[:punct:]]", " ", txtMerged) # remove punctuation with gsub
  txtMerged <- gsub("[^[:alnum:]]", " ", txtMerged) # remove all non-alphanumeric characters

  corpus <- tm::Corpus(tm::VectorSource(txtMerged))
  corpus <- tm::tm_map(corpus, tm::removePunctuation) # remove punctuation with removePunctuation from package tm
  tdm <- tm::TermDocumentMatrix(corpus)
  m <- as.matrix(tdm)
  d <- data.frame(freq = sort(rowSums(m), decreasing = TRUE))

  # Stem words
  d$stem <- SnowballC::wordStem(row.names(d), language = "english")
  d$word <- row.names(d)

  # remove short and long string:
  d <- d[nchar(row.names(d)) < maxword, ]
  d <- d[nchar(row.names(d)) > minword, ]
  d <- d[!grepl(row.names(d), pattern = "www"),] # remove words containing "www"
  d <- d[!grepl(row.names(d), pattern = "http"),] # remove words containing "http"

  # filter with word freq
  d <- d[d$freq >= minFreqWord,]

  if (length(d[,1]) > 0){
    agg_freq <- stats::aggregate(freq ~ stem, data = d, sum) ###
    agg_word <- stats::aggregate(word ~ stem, data = d, function(x) x[1]) ###
    d <- cbind(freq = agg_freq[,2], agg_word) ###
    # sort by freq
    d <- d[order(d$freq, decreasing = T), ]
  }else{}
  return(d)
}
