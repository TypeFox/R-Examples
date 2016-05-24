#' Out-of-Range Report
#'
#' Generate a report on the frequency of variables found outside
#' the user-defined range.
#'
#' @param data data.frame. Contains information required to perform checks.
#' \sQuote{dataframe} should be the name of the data.frame to check.
#' \sQuote{variable} is the name of the variable to check.
#' \sQuote{label} is the variable label.
#' \sQuote{min} is the minimum value.
#' \sQuote{max} is the maximum value.
#' \sQuote{units} is the unit of measurement.
#' @param colheader character. Column header for table.
#' @param panel character. Name for panel.
#' @param append append logical. If \sQuote{TRUE} output will be appended instead of overwritten.
#' @export
#' @examples
#' \dontrun{
#'   load(url('http://biostat.mc.vanderbilt.edu/wiki/pub/Main/Rreport/ssafety.rda'))
#'   rules <- data.frame(
#'     dataframe = rep('ssafety', 4),
#'     variable = c('age', 'height', 'weight', 'bmi'),
#'     label = c('age', 'height', 'weight', 'bmi'),
#'     min = c(45, 145, 50, 15),
#'     max = c(80, 180, 140, 40),
#'     units = c('years', 'cms', 'kgs', 'cm/kg')
#'   )
#'   rangeCheck(rules, panel='check')
#' }

rangeCheck <- function(data, colheader = 'Variable',
                                        panel, append=FALSE) {

  x <- data
  # Make sure valid column names were specified for each dataframe specified in data
  illnames <- NULL
  for(i in as.character(unique(x$dataframe))) {
    if(any(subset(x, dataframe == i)$variable %nin% names(get(i)))) illnames <- c(illnames, i)
  }
  if(length(illnames)) stop(paste('Illegal variable names specified for', paste(illnames, collapse = ', ')))

  # Build a table which specifies the variable, its defined range (including units), and the frequency
  #	of values outside the defined range (% and raw frequency)
  Table <- data.frame(column = NA, min = NA, max = NA, out1 = NA, out2 = NA)
  for(i in 1:nrow(x)) {
    vec <- get(as.character(x[i, 'dataframe']))[ as.character(x[i, 'variable']) ]
    n <- length(vec[!is.na(vec)])
    check <- paste('datta <', x[i, 'min'], '| datta >', x[i, 'max'])
    checkfn <- eval(parse(text = paste('function(datta) {', check, '}')))
    bad <- checkfn(vec)
    nbad <- sum(bad, na.rm = TRUE)
    # Split 'pctbad' by the decimal point so can align the column by the decimal --> see col.just
    pctbad <- unlist(strsplit(format(round((nbad/n)*100, 2), nsmall = 2), '\\.'))
    Table[i, 'column'] <- as.character(x[i, 'label'])
    Table[i, 'min'] <- format(x[i, 'min'], big.mark = ',')
    Table[i, 'max'] <- paste(format(x[i, 'max'], big.mark = ','),
      latexTranslate(as.character(x[i, 'units']), greek = TRUE))
    Table[i, 'out1'] <- pctbad[1]
    Table[i, 'out2'] <- paste(pctbad[2], '\\% \\scriptsize $\\frac{', nbad, '}{', n, '}$', sep = '')
  }

  invisible(latex(Table,
    file = file.path(TexDirName(), sprintf("%s.tex", panel)),
    where = '!htbp', ctable = TRUE, append = FALSE, rowname = NULL,
    cgroup = c(paste('\\textbf{', colheader, '}'), '\\textbf{Defined Range}', '\\textbf{Out of Range}'), 
    n.cgroup = c(1, 2, 2), colheads = NULL,
    # Align the range column by the '-' and the percent bad column by the '.'
    # --> NOTE: col.just must have the same number of elements as ncol(Table) --> have to combine the justifications
    #		to align by the '-' and '.'
    col.just = Cs('l', 'r@{~-~}', 'l', 'r@{.}', 'l'),
    caption = paste('Frequency of', casefold(colheader), 
      'values outside of defined ranges.'))
    )
}
