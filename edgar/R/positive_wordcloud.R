#' Creates wordcloud of positive words from 10-K statement.
#'
#' \code{PositiveWordcloud} creates the wordcloud of positive words from 10-K statement.
#'
#' PositiveWordcloud function takes words frequency dataframe as an input from 
#' \link[edgar]{GetWordfrquency} function. It compares the words frequency dataframe 
#' with the positive words mentioned in the 
#' Loughran and McDonald's financial sentiment dictionaries
#' and generates wordcloud using only positive words with their frequencies.
#'  
#' @param word.frq Word frequency dataframe created using 
#' \link[edgar]{GetWordfrquency} function.
#' 
#' @return Function creates wordcloud of positive words containing in the selected
#' 10-K statement.
#'   
#' @examples
#' \dontrun{
#' 
#' PositiveWordcloud(word.frq)
#' }
#' 
#' @references Loughran and McDonald's financial Sentiment dictionaries
#' \url{http://www3.nd.edu/~mcdonald/Word_Lists.html}
#' @references Bill McDonald, and Tim Loughran.
#' Measuring Readability in Financial Disclosures.
#' Journal of Finance: Volume 69, Issue 4, August 2014
#' @references Bill McDonald, and Tim Loughran.
#' When Is a Liability Not a Liability? Textual Analysis, Dictionaries, and 10-Ks.
#' Journal of Finance: Volume 66, Issue 1, February 2011

PositiveWordcloud <- function(word.frq) {
    
    if (!is.data.frame(word.frq)) {
        msg1 <- "Word_frequency is not a dataframe"
        err <- tcltk::tkmessageBox(message = msg1, icon = "error")
        stop(msg1)
    }
    
    if (nrow(word.frq) == 0) {
        msg2 <- "Word_frequency dataframe is empty"
        err <- tcltk::tkmessageBox(message = msg2, icon = "error")
        stop(msg2)
    }
    
    col <- paste0(names(word.frq), collapse = " ")
    
    if (grepl("FREQUENCY", col) && grepl("WORD", col)) {
		# Read positive words dictionary
        pos.words <- read.csv(system.file("data/poswords.csv", package = "edgar"))
        pos.words <- pos.words$WORDS
        words <- unlist(word.frq$WORD)
        pos.word.table <- word.frq[words %in% pos.words, ]
		
		# Creates wordcloud
        wordcloud::wordcloud(words = pos.word.table$WORD, freq = pos.word.table$FREQUENCY, 
						     scale = c(4, 0.8), max.words = Inf, random.order = F, 
                             colors = RColorBrewer::brewer.pal(8, "Dark2"))
        
    } else {
        msg2 <- "Word_frequency dataframe is invalid"
        err <- tcltk::tkmessageBox(message = msg2, icon = "error")
        stop(msg2)
    }
    
} 
