#' Creates wordcloud of negative words from 10-K statement.
#'
#' \code{NegativeWordcloud} creates the wordcloud of negative words from 10-K statement.
#'
#' NegativeWordcloud function takes words frequency dataframe as an input from 
#' \link[edgar]{GetWordfrquency} function. It compares this words frequency dataframe 
#' with the negative words mentioned in the 
#' Loughran and McDonald's financial sentiment dictionaries
#' and generates wordcloud using only negative words with their frequencies.
#'  
#' @param word.frq Word frequency dataframe created using 
#' \link[edgar]{GetWordfrquency} function.
#' 
#' @return Function creates wordcloud of negative words containing in the selected
#' 10-K statement.
#'   
#' @examples
#' \dontrun{
#' 
#' NegativeWordcloud(word.frq)
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

NegativeWordcloud <- function(word.frq) {
    
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
        # read negative word dictionary
        neg.words <- utils::read.csv(system.file("data/negwords.csv", package = "edgar"))
        neg.words <- neg.words$WORDS
        
        words <- unlist(word.frq$WORD)
        neg.word.table <- word.frq[words %in% neg.words, ]
		# Creates wordcloud
        wordcloud::wordcloud(words = neg.word.table$WORD, freq = neg.word.table$FREQUENCY, 
							 scale = c(4, 0.8), max.words = Inf, random.order = F, 
                             colors = RColorBrewer::brewer.pal(8, "Dark2"))
    } else {
        msg2 <- "Word_frequency dataframe is invalid"
        err <- tcltk::tkmessageBox(message = msg2, icon = "error")
        stop(msg2)
    }
    
} 
