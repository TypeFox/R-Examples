#' @title clean data
#' 
#' @description remove Punctuations, remove Numbers, Translate characters to lower or upper case, remove stopwords, remove user specified words, Stemming words.
#' @param abstracts output of getAbstracts, or just a paragraph of text
#' @param rmNum     Remove the text document with any numbers in it or not
#' @param tolw      Translate characters in character vectors to lower case or not
#' @param toup      Translate characters in character vectors to upper case or not
#' @param rmWords   Remove a set of English stopwords (e.g., 'the') or not
#' @param yrWords   A character vector listing the words to be removed.
#' @param stemDoc   Stem words in a text document using Porter's stemming algorithm.
#' @seealso \code{\link{getAbstracts}}
#' @export
#' @examples
#' # Abs=getAbstracts(c("22693232", "22564732"))
#' # cleanAbs=cleanAbstracts(Abs)
#' 
#' # text="Jobs received a number of honors and public recognition." 
#' # cleanD=cleanAbstracts(text)
cleanAbstracts <- function(abstracts,rmNum=TRUE,tolw=TRUE,toup=FALSE,
                           rmWords=TRUE,yrWords=NULL,stemDoc=FALSE){
  abstTxt <- Corpus(VectorSource(abstracts))
  text2.corpus = tm_map(abstTxt, removePunctuation)
  if(rmNum==TRUE){
    text2.corpus = tm_map(text2.corpus, function(x) removeNumbers(x))
  }
  if(tolw==TRUE){
    text2.corpus = tm_map(text2.corpus, tolower)
  }
  if(toup==TRUE){
    text2.corpus = tm_map(text2.corpus, toupper)
  }
  if(rmWords==TRUE){
    text2.corpus = tm_map(text2.corpus, removeWords, stopwords("english"))
    if(!is.null(yrWords)){
      text2.corpus = tm_map(text2.corpus, removeWords, yrWords)
    }
  }
  if(stemDoc==TRUE){
    text2.corpus = tm_map(text2.corpus, stemDocument)
  }
  text2.corpus <- tm_map(text2.corpus, PlainTextDocument) ### new added
  tdm <- TermDocumentMatrix(text2.corpus)
  m <- as.matrix(tdm)
  v <- sort(rowSums(m),decreasing=TRUE)
  d <- data.frame(word = names(v),freq=v)
  return(d)
}
