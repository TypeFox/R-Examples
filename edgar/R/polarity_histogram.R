#' Creates polarity histogram of word frequency dataframe of 10-K statement.
#'
#' \code{PolarityHist} creates polarity histogram of word 
#' frequency dataframe of 10-K statement.
#'
#' PolarityHist function takes words frequency dataframe as an input from 
#' \link[edgar]{GetWordfrquency} function. It compares this words frequency dataframe 
#' with the negative and positive words mentioned in the 
#' Loughran and McDonald's financial sentiment dictionaries
#' and generates polarity histogram with their frequencies.
#'  
#' @param word.frq Word frequency dataframe created using 
#' \link[edgar]{GetWordfrquency} function.
#' 
#' @return Function creates polarity histogram.
#'   
#' @examples
#' \dontrun{
#' 
#' PolarityHist(word.frq)
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

PolarityHist <- function(word.frq) {
  
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
    
    # read positive dictionary
    pos.words <- utils::read.csv(system.file("data/poswords.csv", package = "edgar"))
    pos.words <- pos.words$WORDS
    pos.word.table <- word.frq[words %in% pos.words, ]
    
    polaritydata <- data.frame(Polarity = c("Negative", "Positive"), 
            Total_words = c(sum(neg.word.table$FREQUENCY), sum(pos.word.table$FREQUENCY))) 
    
    # Creates polarity histogram
    with(polaritydata,{
    ggplot2::ggplot(data = polaritydata, ggplot2::aes(x = Polarity, y = Total_words)) + 
    ggplot2::geom_bar(fill = c("#FF0000", "#339900"), stat = "identity", width = 0.4, position = ggplot2::position_dodge(width = 0.5)) + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(hjust = 0.5)) + 
    ggplot2::theme(axis.text = ggplot2::element_text(size = 14, face = "bold"), axis.title = ggplot2::element_text(size = 16, face = "bold")) + 
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white", colour = "black")) + 
    ggplot2::xlab("Polarity") + 
    ggplot2::ylab("Frequency") + 
    ggplot2::theme(plot.title = ggplot2::element_text(lineheight = 3, face = "bold", color = "black", size = 18))
    })
  } else {
    msg2 <- "Word_frequency dataframe is invalid"
    err <- tcltk::tkmessageBox(message = msg2, icon = "error")
    stop(msg2)
  }
  
}
