#' Read sequences from .txt file
#'
#' Read sequence data saved in text file.
#'
#' @param connection a \code{\link{connection}} to the text (.txt) file.
#' @keywords manip
#' @return a list of sequences. Each element has class \code{\link[seqinr]{SeqFastaAA}}. If
#' connection contains no characters, function prompts warning and returns \code{NULL}.
#' @details The input file should contain one or more amino acid sequences separated by 
#' empty line(s).
#' @export
#' @keywords manip

read_txt <- function(connection) {
  content <- readLines(connection)
  
  #test for empty content
  if(content[1] != "" || length(content) > 1) {
    if (sum(grepl(">", content, fixed = TRUE)) == 0) {
      if (content[1] != "")
        content <- c("", content)
      
      #number of empty lines
      nel <- 0
      #content without too many empty lines
      content2 <- c()
      for (i in 1L:length(content)) {
        if(content[i] == "") {
          nel <- nel + 1
        } else {
          nel <- 0
        }
        if (nel <= 1)
          content2 <- c(content2, content[i])
      }
      content <- content2
      content_end <- length(content)
      while(content[content_end] == "i")
        content_end <- content_end - 1
      prot_names <- sapply(1L:sum(content == ""), function(i)
        paste0(">sequence", i))
      content[content == ""] <- prot_names
    }
    read.fasta(textConnection(content), seqtype = "AA", as.string = FALSE)
  } else {
    warning("No text detected.")
    NULL
  } 
}