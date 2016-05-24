#' Creates histogram of most frequent negative words in 10-K statement.
#'
#' \code{NegWordsHist} creates histogram of most frequent 
#' negative words in 10-K statement.
#'
#' NegWordsHist function takes words frequency dataframe as an input from 
#' \link[edgar]{GetWordfrquency} function. It compares this words 
#' frequency dataframe with negative words mentioned in the 
#' Loughran and McDonald's financial sentiment dictionaries
#' and generates histogram of 15 most frequent negative words 
#' with their frequencies.
#'
#' @param word.frq Word frequency dataframe created using 
#' \link[edgar]{GetWordfrquency} function.
#' 
#' @return Function creates histogram.
#'   
#' @examples
#' \dontrun{
#' 
#' NegWordsHist(word.frq)
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

NegWordsHist <- function(word.frq) {
    
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
        words <- unlist(word.frq$WORD)
        # read negative dictionary
        neg.words <- utils::read.csv(system.file("data/negwords.csv", package = "edgar"))
        neg.words <- neg.words$WORDS
        neg.word.table <- word.frq[words %in% neg.words, ]
        # get top 15 negative words occurred in 10-K statement
        neghistdata <- neg.word.table[1:15, ]
		
		# Creates negative words histogram
        ggplot2::qplot(neghistdata$WORD, weight = neghistdata$FREQUENCY, data = neghistdata, geom = "bar", xlab = "Word", ylab = "Frequency") + 
        ggplot2::coord_flip() + ggplot2::theme(axis.text.x = ggplot2::element_text(hjust = 1, vjust = 0.5)) + ggplot2::geom_bar(fill = "#FF0000") + 
		ggplot2::theme(axis.text = ggplot2::element_text(size = 11, face = "bold"), axis.title = ggplot2::element_text(size = 13, face = "bold")) + 
		ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white", colour = "black")) + ggplot2::ggtitle("Negative words Histogram") + 
		ggplot2::theme(plot.title = ggplot2::element_text(lineheight = 3, face = "bold", color = "black", size = 18))
        
    } else {
        msg2 <- "Word_frequency dataframe is invalid"
        err <- tcltk::tkmessageBox(message = msg2, icon = "error")
        stop(msg2)
    }
    
} 
