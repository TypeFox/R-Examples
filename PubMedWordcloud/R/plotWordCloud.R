#' @title PubMed wordcloud using function 'wordcloud' of package {wordcloud}
#' 
#' @param abs output of cleanAbstracts, or a data frame with one colume of 'word' and one colume of 'freq'.
#' @param scale A vector of length 2 indicating the range of the size of the words.
#' @param min.freq words with frequency below min.freq will not be plotted
#' @param max.words Maximum number of words to be plotted. least frequent terms dropped
#' @param random.order plot words in random order. If false, they will be plotted in decreasing frequency
#' @param rot.per proportion words with 90 degree rotation
#' @param use.r.layout if false, then c++ code is used for collision detection, otherwise R is used
#' @param colors color words from least to most frequent
#' @details This function just call 'wordcloud' from package {wordcloud}. See package {wordcloud} for more details about the parameters.
#' @export
#' @examples
#' # text="Jobs received a number of honors and public recognition." 
#' # cleanD=cleanAbstracts(text)
#' # plotWordCloud(cleanD,min.freq=1,scale=c(2,1))
plotWordCloud <- function(abs, scale=c(3,0.3),min.freq=1,max.words=100,
                          random.order=FALSE,rot.per=0.35,use.r.layout=FALSE, 
                          colors=brewer.pal(8,"Dark2")){
  abs$word=as.character(abs$word)
  abs$freq=as.numeric(as.character(abs$freq))
  abs=abs[abs$freq >= min.freq,]
  wordcloud(abs$word,abs$freq, scale=scale, min.freq = min.freq, max.words=max.words, random.order = random.order, rot.per=rot.per, use.r.layout=use.r.layout,colors=colors)   
}
