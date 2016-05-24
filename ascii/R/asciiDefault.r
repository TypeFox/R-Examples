##' @param x An R object of class found among
##' \code{methods(ascii)}. If \code{x} is a list, it should be a list
##' of character strings (it will produce a bulleted list output by
##' default).
##' @param include.rownames logical. If \code{TRUE} the rows names are printed.
##'   Default value depends of class of \code{x}.
##' @param include.colnames logical. If \code{TRUE} the columns names are
##'   printed. Default value depends of class of \code{x}.
##' @param rownames Character vector (replicated or truncated as necessary)
##'   indicating rownames of the corresponding rows. If \code{NULL} (default)
##'   the row names are not modified
##' @param colnames Character vector (replicated or truncated as necessary)
##'   indicating colnames of the corresponding columns. If \code{NULL}
##'   (default) the column names are not modified
##' @param format Character vector or matrix indicating the format for the
##'   corresponding columns.  These values are passed to the \code{formatC}
##'   function.  Use \code{"d"} (for integers), \code{"f"}, \code{"e"},
##'   \code{"E"}, \code{"g"}, \code{"G"}, \code{"fg"} (for reals), or
##'   \code{"s"} (for strings).  \code{"f"} gives numbers in the usual
##'   \code{xxx.xxx} format; \code{"e"} and \code{"E"} give \code{n.ddde+nn} or
##'   \code{n.dddE+nn} (scientific format); \code{"g"} and \code{"G"} put
##'   \code{x[i]} into scientific format only if it saves space to do so.
##'   \code{"fg"} uses fixed format as \code{"f"}, but \code{digits} as number
##'   of \emph{significant} digits.  Note that this can lead to quite long
##'   result strings. Finaly, \code{"nice"} is like \code{"f"}, but with 0
##'   digits if \code{x} is an integer. Default depends on the class of
##'   \code{x}.
##' @param digits Numeric vector of length equal to the number of columns of
##'   the resulting table (otherwise it will be replicated or truncated as
##'   necessary) indicating the number of digits to display in the
##'   corresponding columns.  Default is \code{2}.
##' @param decimal.mark The character to be used to indicate the numeric
##'   decimal point.  Default is \code{"."}.
##' @param na.print The character string specifying how \code{NA} should be
##'   formatted specially. Default is "".
##' @param caption Character vector of length 1 containing the table's caption
##'   or title.  Set to \code{""} to suppress the caption.  Default value is
##'   \code{NULL}.
##' @param caption.level Character or numeric vector of length 1 containing the
##'   caption's level.  Can take the following values: \code{0} to \code{5},
##'   \code{"."} (block titles in asciidoc markup), \code{"s"} (strong),
##'   \code{"e"} (emphasis), \code{"m"} (monospaced) or \code{""} (no markup).
##'   Default is NULL.
##' @param width Numeric vector of length one containing the table width
##'   relative to the available width (expressed as a percentage value,
##'   \code{1}\dots{} \code{99}).  Default is \code{0} (all available width).
##' @param frame Character vector of length one. Defines the table border, and
##'   can take the following values: \code{"topbot"} (top and bottom),
##'   \code{"all"} (all sides), \code{"none"} and \code{"sides"} (left and
##'   right).  The default value is \code{NULL}.
##' @param grid Character vector of length one. Defines which ruler lines are
##'   drawn between table rows and columns, and can take the following values:
##'   \code{"all"}, \code{"rows"}, \code{"cols"} and \code{"none"}.  Default is
##'   \code{NULL}.
##' @param valign Vector or matrix indicating vertical alignment of all cells
##'   in table.  Can take the following values: \code{"top"}, \code{"bottom"}
##'   and \code{"middle"}.  Default is \code{""}.
##' @param header logical or numeric. If \code{TRUE} or \code{1}, \code{2},
##'   \dots{}, the first line(s) of the table is (are) emphasized. The default
##'   value depends of class of \code{x}.
##' @param footer logical or numeric. If \code{TRUE} or \code{1}, the last
##'   line(s) of the table is (are) emphasized. The default value depends of
##'   class of \code{x}.
##' @param align Vector or matrix indicating the alignment of the corresponding
##'   columns.  Can be composed with \code{"r"} (right), \code{"l"} (left) and
##'   \code{"c"} (center).  Default value is \code{NULL}.
##' @param col.width Numeric vector of length equal to the number of columns of
##'   the resulting table (otherwise it will be replicated or truncated as
##'   necessary) indicating width of the corresponding columns (integer
##'   proportional values).  Default is \code{1}.
##' @param style Character vector or matrix indicating the style of the
##'   corresponding columns.  Can be composed with \code{"d"} (default),
##'   \code{"s"} (strong), \code{"e"} (emphasis), \code{"m"} (monospaced),
##'   \code{"h"} (header) \code{"a"} (cells can contain any of the AsciiDoc
##'   elements that are allowed inside document), \code{"l"} (literal),
##'   \code{"v"} (verse; all line breaks are retained).  Default is
##'   \code{NULL}.
##' @param tgroup Character vector or a list of character vectors defining
##'   major top column headings.  The default is to have none (\code{NULL}).
##' @param n.tgroup A numeric vector or a list of numeric vectors containing
##'   the number of columns for which each element in tgroup is a heading.  For
##'   example, specify \code{tgroup=c("Major 1","Major 2")},
##'   \code{n.tgroup=c(3,3)} if \code{"Major 1"} is to span columns 1-3 and
##'   \code{"Major 2"} is to span columns 4-6.
##' @param talign Character vector of length one defining alignment of major
##'   top column headings.
##' @param tvalign Character vector of length one defining vertical alignment
##'   of major top column headings.
##' @param tstyle Character vector of length one indicating the style of major
##'   top column headings
##' @param bgroup Character vector or list of character vectors defining major
##'   bottom column headings.  The default is to have none (\code{NULL}).
##' @param n.bgroup A numeric vector containing the number of columns for which
##'   each element in bgroup is a heading.
##' @param balign Character vector of length one defining alignment of major
##'   bottom column headings.
##' @param bvalign Character vector of length one defining vertical alignment
##'   of major bottom column headings.
##' @param bstyle Character vector of length one indicating the style of major
##'   bottom column headings
##' @param lgroup Character vector or list of character vectors defining major
##'   left row headings.  The default is to have none (\code{NULL}).
##' @param n.lgroup A numeric vector containing the number of rows for which
##'   each element in lgroup is a heading. Column names count in the row
##'   numbers if \code{include.colnames = TRUE}.
##' @param lalign Character vector of length one defining alignment of major
##'   left row headings.
##' @param lvalign Character vector of length one defining vertical alignment
##'   of major left row headings.
##' @param lstyle Character vector of length one indicating the style of major
##'   left row headings
##' @param rgroup Character vector or list of character vectors defining major
##'   right row headings.  The default is to have none (\code{NULL}).
##' @param n.rgroup A numeric vector containing the number of rows for which
##'   each element in rgroup is a heading. Column names count in the row
##'   numbers if \code{include.colnames = TRUE}.
##' @param ralign Character vector of length one defining alignment of major
##'   right row headings.
##' @param rvalign Character vector of length one defining vertical alignment
##'   of major right row headings.
##' @param rstyle Character vector of length one indicating the style of major
##'   right row headings
##' @param list.type Character vector of length one indicating the list type
##'   (\code{"bullet"}, \code{"number"}, \code{"label"} or \code{"none"}). If
##'   \code{"label"}, \code{names(list)} is used for labels. Default is
##'   \code{"bullet"}.
##' @param ... Additional arguments. (Currently ignored.)
##' @keywords print
##' @rdname ascii
##' @export
##' @method ascii default
ascii.default <- function(x, include.rownames = TRUE, include.colnames = TRUE, rownames = NULL, colnames = NULL, format = "f", digits = 2, decimal.mark = ".", na.print = "", caption = NULL, caption.level = NULL, width = 0, frame = NULL, grid = NULL, valign = NULL, header = TRUE, footer = FALSE, align = NULL, col.width = 1, style = NULL, tgroup = NULL, n.tgroup = NULL, talign = "c", tvalign = "middle", tstyle = "h", bgroup = NULL, n.bgroup = NULL, balign = "c", bvalign = "middle", bstyle = "h", lgroup = NULL, n.lgroup = NULL, lalign = "c", lvalign = "middle", lstyle = "h", rgroup = NULL, n.rgroup = NULL, ralign = "c", rvalign = "middle", rstyle = "h", list.type = "bullet", ...) {
  if (is.list(x)) {
    x <- lapply(x, as.character)
    obj <- asciiList$new(x = x, caption = caption, caption.level = caption.level, list.type = list.type)
  }
  else {
    y <- as.data.frame(x)
    obj <- asciiTable$new(x = y, include.rownames = include.rownames, include.colnames = include.colnames, rownames = rownames, colnames = colnames, format = format, digits = digits, decimal.mark = decimal.mark, na.print = na.print, caption = caption, caption.level = caption.level, width = width, frame = frame, grid = grid, valign = valign, header = header, footer = footer, align = align, col.width = col.width, style = style, tgroup = tgroup, n.tgroup = n.tgroup, talign = talign, tvalign = tvalign, tstyle = tstyle, bgroup = bgroup, n.bgroup = n.bgroup, balign = balign, bvalign = bvalign, bstyle = bstyle, lgroup = lgroup, n.lgroup = n.lgroup, lalign = lalign, lvalign = lvalign, lstyle = lstyle, rgroup = rgroup, n.rgroup = n.rgroup, ralign = ralign, rvalign = rvalign, rstyle = rstyle)
  }

  return(obj)
}
