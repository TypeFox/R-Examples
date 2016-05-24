#' Export the transformation equation into LaTeX
#' 
#' \code{transeq} turns a \code{jtrans} object into a LaTeX equation for 
#' display.
#' 
#' A LaTeX equation in the display mode, e.g. between \code{\\[} and \code{\\]}
#' is returned with the formula used in the transformation. Note that when it's
#' displayed in the R console, the backslashes are escaped. So it's always
#' double backslash when in print it in the terminal.
#' 
#' This is designed to work with \pkg{knitr} and \pkg{rmarkdown}. In this
#' case you can set the chunk option \code{results='asis'} and output it to a
#' PDF document. Then the LaTeX equation will be properly formatted and can be
#' easily included in your report.
#' 
#' @param obj a \code{jtrans} object with a specific transformation type
#' @param digits digits displayed in the equation
#'
#' @export
#' @examples 
#' 
#' \dontrun{
#' # designed to be used with R Markdown and chunk options
#' ```{r, results='asis'}
#' library(jtrans)
#' jt <- jtrans(rexp(30, .3))
#' transeq(jt)
#' ```
#' }
#' @rdname transeq
transeq <- function(obj, digits = 4) {
  UseMethod("transeq", obj)
}

#' @rdname transeq
#' @export
transeq.sb <- function(obj, digits) {
  form <- paste("\\[ Y=\\gamma+\\eta \\ln\\frac{X-\\epsilon}{", 
                "\\lambda+\\epsilon-x} = ",
                format(obj$gamma, digits = digits), "+(", 
                format(obj$eta, digits = digits), ")\\ln\\frac{X-(", 
                format(obj$epsilon, digits = digits), 
                ")}{", format(obj$lambda, digits = digits), "+(", 
                format(obj$epsilon, digits = digits), ")-X} \\]")
  cat(form)
}

#' @export
transeq.sl <- function(obj, digits) {
  form <- paste("\\[ Y = \\gamma+\\eta \\ln(X-\\epsilon) = ",
                format(obj$gamma, digits = digits), "+(", 
                format(obj$eta, digits = digits), ")\\ln(X-(", 
                format(obj$epsilon, digits = digits), ") \\]")
  cat(form)
}

#' @rdname transeq
#' @export
transeq.su <- function(obj, digits) {
  form <- paste("\\[ Y=\\gamma+\\eta", 
                "\\textrm{sinh}^{-1}\\frac{x-\\epsilon}{\\lambda}=",
                format(obj$gamma, digits = digits), "+(", 
                format(obj$eta, digits = digits), 
                ")\\textrm{sinh}^{-1}\\frac{X-(", 
                format(obj$epsilon, digits = digits), 
                ")}{", format(obj$lambda, digits = digits), "} \\]")
  cat(form)
}



