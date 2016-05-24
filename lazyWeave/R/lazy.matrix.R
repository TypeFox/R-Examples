#' @name lazy.matrix
#' @export lazy.matrix
#' 
#' @title Convert Matrix to LaTeX Table
#' @description An example of using \code{lazyWeave} to produce tables
#' 
#' @param x A matrix.  Other objects are coerced to matrices
#' @param align Character vector or string giving the alignment for each 
#'   column of the table.  Options are \code{"left", "center", "right"}.
#' @param justify Character string giving the alignment for the table on the 
#'   page. Options are \code{"left", "center", "right"}.
#' @param rcol A vector giving the rows of \code{x} to be colored
#' @param usecol A character string or vector denoting the color of the rows
#'   in \code{rcol}
#' @param caption Caption for the table.  This is printed above the table
#' @param footnote Additional footnotes for the table.  These are printed
#'   below the table.
#' @param placement Controls the placement of the figure.  Options are
#'   \code{"ht", "t", "b", "p", "H"} and can be supplemented with 
#'   \code{"!"}. See "Details" for more explanation. These apply only to LaTeX
#' @param translate Toggles if inputs in \code{x} should be passed through 
#'   \code{latexTranslate}.  This should be set to \code{FALSE} if writing
#'   custom code
#' @param cat Logical. Determines if the output is returned as a character string
#'   or returned via the \code{cat} function (printed to console).  The default
#'   value is set by \code{options()$lazyWeave_cat}.  This argument allows for
#'   selective override of the default.
#' @param ... Additional arguments to be passed to \code{lazy.table}
#' 
#' @details The output for \code{lazy.matrix} is highly inflexible compared 
#' to \code{lazy.table}.  It is an example of how to build a reproducible 
#' table with a certain formatting style that may be used again and again
#' for consistency.
#' 
#' Row names are always left justified.  This cannot be changed.
#' 
#' \code{placement} options are used as follows:
#' \tabular{ll}{
#'     ht \tab Place the float here, i.e., 
#'     approximately at the same point it occurs \cr
#'     t  \tab Position at the top of the page\cr
#' b  \tab Position at the bottom of the page \cr
#'     p  \tab Put on a special page for floats only \cr
#'     H  \tab Places the float at precisely the location in the LaTeX code. 
#'     Requires the float package\cr
#'   }
#' The \code{"!"} may be used after any of these in order to override 
#' LaTeX float rules and force your selection.  More can be learned by 
#' reading about floats in a LaTeX manual
#' 
#' @author Benjamin Nutter
#' 
#' @examples 
#' \dontrun{
#' lazy.write(
#'   lazy.file.start(),
#'   lazy.text("The mtcars dataset describes a number of vehicles.  
#'       Let's take a look at the data"),
#'   lazy.matrix(mtcars, rcol=(1:nrow(mtcars))[c(FALSE, TRUE)]),
#'   lazy.file.end(),
#'   OutFile="Example 1.tex")
#' 
#' unlink("Example 1.tex")
#'}
#'


lazy.matrix <-
function(x, align="center", justify="center", rcol=NULL, usecol="lightgray",
    caption=NULL, footnote=NULL, placement="h", translate=TRUE, 
    cat=getOption("lazyWeave_cat"), ...){
    
     
  #*** retrieve the report format
  reportFormat <- getOption("lazyReportFormat")
  if (!reportFormat %in% c("latex", "html", "markdown")) stop("option(\"lazyReportFormat\") must be either 'latex', 'html', or 'markdown'")
  
  #*** Construct the comment with the function call
  comment.char <- if (reportFormat == "latex") {
    if (getOption("lazyWeave_latexComments") == "latex") c("%%", "") else c("<!-- ", " -->")
  }  else if (reportFormat == "html") c("<!--", "-->")
  
  fncall <- paste(comment.char[1], paste(deparse(match.call()), collapse=" "), comment.char[2], "\n")

#*** Coerce x to a matrix
  if (!is.matrix(x)) x <- as.matrix(x)
  
  if ("cwidth" %in% names(list(...))){
    cw <- list(...)$cwidth
    if (!is.null(rownames(x))){
      if (length(cw) != 1 && ((ncol(x) + 1) != length(cw)))
        stop("'cwidth' must have length 1 or equal to ncol(x)--remember your row names")
    }
  }
  
#*** Extend length of align to number of columns of x.  This will be useful
#*** if we add a column for rownames
  if (length(align) == 1) align <- rep(align, ncol(x))

#*** Add the rownames to x and assign them left justification
  if (!is.null(rownames(x))){
    x <- cbind(rownames(x), x)
    rownames(x) <- NULL
    align = c("left", align)
  }
 

#*** Table if colnames are present
  if (!is.null(colnames(x))){
    colnames(x)[-1] <- lazy.text.format(colnames(x)[-1], bold=TRUE)
    header <- lazy.table(colnames(x), align=align, cspan=1,
                          justify=justify, rborder=c(0, 0, 1), 
                          open=TRUE, close=FALSE,
                          caption=caption, placement=placement,
                          translate=translate, cat=FALSE, ...)
    body <- lazy.table(x, align=align, cspan=1,
                        rborder=nrow(x), rcol=rcol,
                        justify=justify, usecol=usecol, 
                        open=FALSE, close=TRUE, 
                        footnote=footnote,
                        translate=translate, cat=FALSE, ...)
  }
  
#*** Table if colnames are not present
  else{
    header <- ""
    body <- lazy.table(x, align=align, cspan=1,
                        justify=justify, rborder=c(0, nrow(x)),
                        open=TRUE, close=TRUE, 
                        rcol=rcol, usecol=usecol,
                        caption=caption, footnote=footnote,
                        placement=placement,
                        translate=translate, cat=FALSE, ...)
  }
  
  code <- paste(if (reportFormat != "markdown") fncall else "", header, body, sep="")

  if (cat) cat(code)
  else return(code)
}

