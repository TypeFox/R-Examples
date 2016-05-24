#' @name textable
#' @author Sven E. Templer
#' @title Table to Latex
#' @description 
#' This function enhances \code{xtable}: It can write the latex code of the
#' table directly to a file and optionally adds a header/footer for the
#' document structure. Also a system command can be given to convert the
#' tex file to a pdf document, for example.
#' @details
#' Example for a system call:\cr
#' \code{cmd = "pdflatex -output-directory /path/to/files/"}
#' @param d Object (will be coerced to data.frame) to transform to a latex table.
#' @param file Character string with output file name. If missing or \code{""},
#' the output is printed to the screen.
#' @param caption Character vector with title of table.
#' @param label Character vector with the latex label or HTML anchor.
#' @param align Character vector with \code{'l'}, \code{'c'}, \code{'r'}
#' for aligning the columns left, centered or right. Length is either one
#' or 1 (for rownames column) + number of columns in \code{d} (even if
#' \code{rownames = FALSE})
#' @param rownames Logical, include row names of \code{d}.
#' @param topcapt Logical, put caption and label before 'tabular'.
#' @param digits Number of digits to print from numeric values.
#' @param as.document Logical. \code{TRUE} to add the document definition to
#' the output. The document class is an article and the package a4paper is
#' included.
#' @param landscape Logical, use a landscape format for wider tables.
#' Only with \code{as.document=TRUE}.
#' @param pt.size Integer from 10 to 13 for the size of the characters.
#' Only with \code{as.document=TRUE}.
#' @param margin Margin between table and page border in cm.
#' Only with \code{as.document=TRUE}.
#' @param cmd A character vector with the system command to apply
#' on the output file. Only if \code{file} is given and \code{as.document}
#' is \code{TRUE}. \code{NULL} or an empty string \link{system} is not called.
#' @param ... Forwarded arguments to \link[xtable]{print.xtable}.
#' @return
#' Returns a character vector invisible. If \code{file} is set, then the
#' content is written to a file. Else it is printed to the console.
#' @seealso
#' \link{xtable} 
#' @examples
#' #
#' 
#' \dontrun{
#' d <- head(trees)
#' dc <- 'R "trees" dataset.'
#' textable(d, rownames=TRUE, digits=4, caption=dc)
#' textable(d, '/tmp/trees.tex', caption=dc, as.document=TRUE, 
#'   cmd='pdflatex --output-directory /tmp')
#' }
#' 
#' #

#' @rdname textable
#' @export textable
textable <- function (
  d, file, caption = NULL, label = NULL, align = NULL, rownames = FALSE, topcapt = TRUE,
  digits = NULL, as.document = FALSE, landscape = FALSE, margin = 2, pt.size = 10, cmd = NULL,
  ...)
{
  
  # replicate align
	if (!is.null(align) && length(align) == 1)
		align <- rep(align, ncol(d)+1)
  
  # get the table
	tex.tab <- capture.output(print(xtable(
    	d, digits=digits, align=align, label=label, caption=caption),
    	include.rownames = rownames, ...))
  
  # store/drop comments
	tex.cmt <- paste(
	  "% output by function 'textable' from package miscset", 
	  as.character(packageVersion("miscset")))
	i.cmt <- grepl("^%", tex.tab)
  tex.cmt <- c(tex.cmt, tex.tab[i.cmt])
  tex.tab <- tex.tab[!i.cmt]
  
  # switch capture position
  tex.cpt <- character()
  if (topcapt) {
    i.cpt <- which(grepl("\\caption", tex.tab))
    if (length(i.cpt)) {
      i.cptn <- which(grepl("}$", tex.tab[-seq(i.cpt)]))[1]
      i.cpt <- seq(i.cpt, length.out = i.cptn)
      tex.cpt <- tex.tab[i.cpt]
      tex.tab <- tex.tab[-i.cpt]
    }
  }
  
  # tex document header
  if (!is.null(cmd) && !as.document) {
    as.document <- TRUE
    message("Set 'as.document' to TRUE because 'cmd' was provided.")
  }
  tex.dochead <- tex.doctail <- character()
  if (as.document) {
    if (!pt.size %in% 10:12)
      stop('pt.size must be 10, 11 or 12.')
    orientation <- if (landscape) 'landscape,' else ''
    tex.dochead <- paste0(
      '\\documentclass[a4paper,', pt.size,
      'pt]{article}\n\\usepackage[a4paper,',
      orientation,
      'margin=',
      margin,
      'cm]{geometry}\n\\begin{document}\n')
    #\\small\n
    tex.doctail <- '\n\\end{document}\n'
  }
  
  # merge
  tex <- c(
    tex.cmt, 
    '',
    tex.dochead,
    tex.tab[1:2],
    if (topcapt) tex.cpt,
    tex.tab[-(1:2)],
    tex.doctail)
  
  # print
	if (missing(file))
		file <- ""
	cat(tex, file=file, sep='\n')
  
  # tex to pdf
  if (as.document && !missing(file) && nchar(file)>0 && !is.null(cmd) && nchar(cmd)>0)
    system(paste(cmd, file), wait = FALSE)
  
  # return
	invisible(tex)
  
}

