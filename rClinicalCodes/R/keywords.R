#' Extract keywords from clinical code lists
#' 
#' This function takes a dataframe of clinical codes and gives a vector of keywords, sorted by frequency
#' 
#' All terms are converted to lower case.
#' Common stopwords, whitespace and punctuation are removed
#' Optional extra_stopwords vector
#' @export
#' @param codelist a dataframe of clinical codes, such as downloaded using the get_ClinicalCodes function
#' @param keyword_column The column of the dataframe to extract keywords from
#' @param extra_stopwords an optional character vector of further stopwords to remove
#' @return character vector of keywords, sorted by frequency in the list
#' @examples \dontrun{
#' # Get codelist from url:
#' angina_codes <- get_ClinicalCodes(
#' url = "https://clinicalcodes.rss.mhs.man.ac.uk/medcodes/article/6/codelist/angina/download/")
#' codelist_keywords(angina_codes, extra_stopwords = c("good", "poor", "[x]"))
#' } 
codelist_keywords <- function(codelist, keyword_column = "description", extra_stopwords = NULL){
    descriptions <- as.character(codelist[[keyword_column]])
    des <- Corpus(VectorSource(descriptions))
    des <- tm_map(des, stripWhitespace)
    des <- tm_map(des, tolower)
    des <- tm_map(des, removeWords, stopwords("english"))
    if(!is.null(extra_stopwords)) des <- tm_map(des, removeWords, extra_stopwords)
    dtm <- DocumentTermMatrix(des)
    findFreqTerms(dtm)
}
