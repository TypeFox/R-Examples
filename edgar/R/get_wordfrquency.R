#' Creates words frequency dataframe of 10-K statement.
#'
#' \code{GetWordfrquency} creates word frequency dataframe of 10-K statement.
#'
#' GetWordfrquency function asks the user to locate 10-K statement which can be download using 
#' \link[edgar]{Downloadfilings} function. Function cleans text of 10-K statement
#' and creates words frequency dataframe. This words frequency dataframe is used
#' in the functions \link[edgar]{PositiveWordcloud},\link[edgar]{PosWordsHist},
#' \link[edgar]{NegativeWordcloud},\link[edgar]{NegWordsHist}, and \link[edgar]{PolarityHist}.
#'  
#' @return Function returns words frequency dataframe.
#'   
#' @examples
#' \dontrun{
#' 
#' word.frq <- GetWordfrquency()
#' }

GetWordfrquency <- function() {
    filepath <- jchoose.files(default = getwd(), caption = "Select 10-K file", multi = FALSE)
    word.frq <- data.frame()
    
    if (grepl("10-K", filepath) && grepl(".txt", filepath)) {
        text <- readLines(filepath)
        text <- paste(text, collapse = " ")
        
        # Extract text from html file
        doc <- XML::htmlParse(text, asText = TRUE)
        text <- XML::xpathSApply(doc, "//text()[not(ancestor::script)][not(ancestor::style)][not(ancestor::noscript)][not(ancestor::form)]", 
                                 XML::xmlValue)
        text <- paste(text, collapse = " ")
        
        # convert into corpus
        text <- tm::Corpus(tm::VectorSource(text))
        # clean text
        cleantext <- function(data.text.corpus) {
            data.text.corpus <- tm::tm_map(data.text.corpus, tm::removePunctuation)  # Remove punctuation marks
            data.text.corpus <- tm::tm_map(data.text.corpus, tm::removeNumbers)  # Remove Numbers
            data.text.corpus <- tm::tm_map(data.text.corpus, tm::stripWhitespace)  # Remove punctuation marks
            data.text.corpus <- tm::tm_map(data.text.corpus, function(x) tm::removeWords(x, tm::stopwords()))  # Remove stop words
            data.text.corpus <- tm::tm_map(data.text.corpus, tm::content_transformer(tolower))  # Convert text to lower case
            return(data.text.corpus)
        }
        text <- cleantext(text)
        word.frq <- tm::termFreq(text[[1]])
        word.frq <- data.frame(WORD = names(word.frq), FREQUENCY = word.frq, row.names = NULL)
		
		# Order dataframe descending on frequency
        word.frq <- word.frq[order(-word.frq$FREQUENCY), ]
        rownames(word.frq) <- NULL
    } else {
        msg3 <- "Please select 10-K file only.."
        err <- tcltk::tkmessageBox(message = msg3, icon = "error")
    }
    return(word.frq)
} 
