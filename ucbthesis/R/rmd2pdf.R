#' @title Render an R Markdown file into a PDF
#' 
#' @description Use knitr and pandoc to convert an R Markdown file into a
#' TeX file, and run xelatex (and biber) on the converted file to produce a PDF.
#'  
#' @param file the location and name of the R Markdown file to be rendered.
#' @param template the location and name of the pandoc latex template to use
#' during the conversion from R Markdown to TeX. 
#' @param biber logical flag indicating if biber (or biblatex) backend should
#' be run after xelatex is called.
#' @param saveTmpFiles logical flag indicating if intermediary files should be 
#' kept after PDF file is created. If \code{FALSE} the files are deleted.
#' 
#' @details 
#' There's no markdown (yet) that allows for LaTeX preamble to be specified
#' inside a .md file. To get around this, we have to redefine a Pandoc latex
#' template using the preamble we'd normally use if we were using the LaTeX or
#' knitr templates. The template, \code{inst/rmarkdown/thesis_template.latex},
#' is a modification of the default template found via
#' \code{pandoc -D latex}. If you need to add more packages to your preamble,
#' e.g. \code{\\usepackage{amsmath}}, modify \code{thesis_template.latex} 
#' accordingly.
#' 
#' Temporary files (e.g. .md's, .log's, .aux's, etc.) are stored in a temporary 
#' (sub)directory, \code{tmp/}.
#' 
#' @export
#' 
#' @return The name of the xelatex rendered PDF.
#' 
#' @examples
#' \dontrun{
#' setwd("inst/rmarkdown")
#' rmd2pdf()
#' }
rmd2pdf <- function(file='thesis.Rmd', template='thesis_template.latex',
                    biber=TRUE, saveTmpFiles=FALSE) {
  
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
  
  if (str_length(Sys.which('pandoc')) == 0) {
    stop(str_c('Must have pandoc installed and accesible from the command line',
               ' to run this function.')
         )
  }
  
  cat("Working in", getwd(), "...\n")
  
  filename <- tail(
                str_split(
                  string = str_split(file, ".(R|r)md", n=2)[[1]][1],
                  pattern = "/")[[1]],
                n=1)
  
  cat("Making ./tmp/ directory to hold temporary files...\n")
  system("mkdir tmp")
  
  mdFile <- str_c("./tmp/", filename, ".md")
  
  cat("knit'ing", file, "into", mdFile, "...\n") 
  knitr::knit(input = file, output = mdFile)
  
  texFile <- str_c("./tmp/", filename, ".tex")
  
  # make pandoc command string
  pandocCmd <- str_c("pandoc ", shQuote(mdFile),
                     " -t latex -o ", shQuote(texFile),
                     " --template=", shQuote(template))
  
  cat("pandoc'ing", mdFile, "into", texFile, "...\n")
  system(pandocCmd)

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
  
  saveExtensions <- c("Rmd", "rmd", "bib", "latex", "cls", "sty", "pdf")
  
  if (xelatexFail) {
    stop("\n\nxelatex failed to compile a PDF from the rendered .tex file...\n")
    stop("Cleaning up temporary files...\n\n")
    cleanUp(saveExtensions)
  } 
  
  if (!saveTmpFiles) cleanUp(saveExtensions)
  
  return(str_c(filename, ".pdf"))
}
