#' @title Render an Rnw file into a PDF
#' 
#' @description Use knitr to convert an Rnw file into xelatex (and biber) 
#' rendered PDF.
#'  
#' @param file the location and name of the Rnw file to be rendered.
#' @param biber logical flag indicating if biber (or biblatex) backend should
#' be run after xelatex is called.
#' @param saveTmpFiles logical flag indicating if intermediary files should be 
#' kept after PDF file is created. If \code{FALSE} the files are deleted.
#' 
#' @details
#' This is just a sequence of system calls that runs knitr to turn the .Rnw into
#' a .tex, and then calls xelatex (and biber) a few more times to make sure 
#' citations and cross-references are correct.
#' 
#' Temporary files (e.g. .tex's, .log's, .aux's, etc.) are stored in a temporary 
#' (sub)directory, \code{tmp/}.
#' 
#' @export
#' 
#' @return The name of the xelatex rendered PDF.
#' 
#' @examples
#' \dontrun{
#' setwd("inst/knitr")
#' rnw2pdf()
#' }
rnw2pdf <- function(file='thesis.Rnw', biber=TRUE, saveTmpFiles=FALSE) {
  
  if (str_length(Sys.which('xelatex')) == 0) {
    stop(str_c('Must have xelatex installed and accesible from the command line',
               ' to run this function.')
    )
  }
  
  if (biber && str_length(Sys.which('biber')) == 0) {
    stop(str_c('Must have biber installed and accesible from the command line',
               ' to run this function.')
    )
  }
  
  cat("Working in", getwd(), "...\n")
  
  filename <- tail(
                str_split(
                  string = str_split(file, ".(R|r)nw", n=2)[[1]][1],
                  pattern = "/")[[1]],
                n=1)
  
  cat("Making ./tmp/ directory to hold temporary files...\n")
  system("mkdir tmp")
  
  texFile <- str_c("./tmp/", filename, ".tex")
  
  cat("knit'ing", file, "into", texFile, "...\n") 
  knitr::knit(input = file, output = texFile)

  # run xelatex 
  xelatexCmd <- str_c("xelatex ", texFile)
  
  cat("xelatex'ing", texFile, "...\n")
  xelatexFail <- system(xelatexCmd)
  
  if (biber) {
    cat("\n\nRunning biber on intermediate files...\n\n")
    system(str_c("biber ", filename))
    
    cat("\n\nRe-running xelatex to clean up cross-references...\n\n")
    xelatexFail <- system(xelatexCmd)
  }
  
  cat("\n\nMove temporary files created from xelatex:\n")
  moveLatexFiles()
  
  saveExtensions <- c("Rnw", "rnw", "bib", "latex", "cls", "sty", "pdf")
  
  if (xelatexFail) {
    stop("\n\nxelatex failed to compile a PDF from the rendered .tex file...\n")
    stop("Cleaning up temporary files...\n\n")
    cleanUp(saveExtensions)
  } 
  
  if (!saveTmpFiles) cleanUp(saveExtensions)
  
  return(str_c(filename, ".pdf"))
}
