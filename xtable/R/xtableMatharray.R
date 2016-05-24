### xtableMatharray object
### To deal with numeric arrays such as a variance-covariance matrix
### From a request by James Curran, 16 October 2015
xtableMatharray <- function(x, caption = NULL, label = NULL,
                            align = NULL, digits = NULL,
                            display = NULL, auto = FALSE,
                            ...) {
  class(x) <- c("xtableMatharray","matrix")
  xtbl <- xtable.matrix(x,
                        caption = caption, label = label, align = align,
                        digits = digits, display = display, auto = auto,
                        ...)
  class(xtbl) <- c("xtableMatharray","xtable","data.frame")
  return(xtbl)
}

print.xtableMatharray <- function(x,
           print.results = TRUE,
           format.args = getOption("xtable.format.args", NULL),
           scalebox = getOption("xtable.scalebox", NULL),
           comment = FALSE,
           timestamp = NULL,
           ...)
{
  class(x) <- c("xtableMatharray","data.frame")
  print.xtable(x, floating = FALSE,
               tabular.environment = 'array',
               include.rownames = FALSE, include.colnames = FALSE,
               hline.after = NULL,
               print.results = print.results,
               format.args = format.args,
               scalebox = scalebox,
               comment = comment,
               timestamp = timestamp,
               ...)
}
