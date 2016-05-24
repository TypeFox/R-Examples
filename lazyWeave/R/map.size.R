#' @name map.size
#' 
#' @title Convert text size specifications between LaTeX and HTML formats
#' @description The default text size specifications in LaTeX are character 
#'     strings, but in HTML, the font size is determined by an integer.  
#'     In order to be able to switch formats in existing code, there needs 
#'     to be some way to map HTML text size to LaTeX text size, and vice versa.  
#'     This is managed by \code{map.size} in the \code{lazyWeave} functions 
#'     so the user won't have to worry about it.
#' 
#' @param x The text size specification to be mapped to the desired format
#' @param reportFormat The format to which \code{x} should be mapped
#' 
#' @details 
#'   Text sizes are mapped according to the following:
#'   \tabular{ll}{
#'     LaTeX        \tab HTML\cr
#'     tiny         \tab 5\cr
#'     scriptsize   \tab 7\cr
#'     footnotesize \tab 8\cr
#'     small        \tab 9\cr
#'     normalsize   \tab 10\cr
#'     large        \tab 12\cr
#'     Large        \tab 14\cr
#'     LARGE        \tab 18\cr
#'     huge         \tab 20\cr
#'     Huge         \tab 24\cr
#' }
#'     
#' @author Benjamin Nutter
#' 
  
map.size <- function(x, reportFormat=getOption("lazyReportFormat")){
  
  size.ref <- data.frame(latex=c("tiny", "scriptsize", "footnotesize", "small", "normalsize", "large",
                                 "Large", "LARGE", "huge", "Huge"),
                         html=c(5, 7, 8, 9, 10, 12, 14, 18, 20, 24),
                         stringsAsFactors=FALSE)
  
  if (reportFormat == "latex"){
    if (gsub("[[:punct:]]", "", x) %in% size.ref$latex) return(x)
    if (is.numeric(x) && x >= 5) html.size <- utils::tail(size.ref$html[size.ref$html <= x], 1) else html.size <- 5
    x <- paste("\\", size.ref$latex[size.ref$html == html.size], sep="")
  }
  
  if (reportFormat == "html"){
    if (is.numeric(x)) return(x)
    if (gsub("[[:punct:]]", "", x) %in% size.ref$latex) x <- size.ref$html[size.ref$latex == gsub("[[:punct:]]", "", x)]
    else{
      x <- get(".HTML.FONT.SIZE.", envir=.GlobalEnv)
      warning("Could not map character size description to html font size.  The option '.HTML.FONT.SIZE.' is used instead")
    }
  }
      
  return(x)
}      

