#' Subject List
#'
#' Generate a LaTeX table from a dataset.
#'
#' @param data data.frame. Data used for report.
#' @param panel character. Name of panel.
#' @param caption character. See \code{\link[Hmisc]{latex}}.
#' @param vname character. Specifies how to generate column headings,
#' either through variable \sQuote{labels} or \sQuote{names}.
#' @param colheads character vector. Column headings for each variable.
#' @param size character. Set LaTeX table font size, see \code{\link[Hmisc]{latex}}.
#' @param longtable logical. See \code{\link[Hmisc]{latex}}.
#' @param landscape logical. See \code{\link[Hmisc]{latex}}.
#' @export
#' @examples
#' \dontrun{
#'   load(url('http://biostat.mc.vanderbilt.edu/wiki/pub/Main/Rreport/ssafety.rda'))
#'   subjectList(ssafety[1:10,1:10], "datalist", vname='names')
#' }

subjectList <- function(data, panel, caption=NULL,
                        vname=c('labels','names'),
                        colheads=NULL,
                        size='smaller',
                        longtable=TRUE, landscape=TRUE) {

  vname <- match.arg(vname)
  if(length(colheads)) lab <- colheads else {
    lab <- names(data)
    if(vname == 'labels') {
      lab <- sapply(data,label)
      lab <- ifelse(lab=='', names(data), lab)
    }
  }
  ## For chron date-time variables remove surrounding ( ) and seconds
  for(i in 1:length(data)) {
    x <- data[[i]]
    if(all(c('chron','dates','times') %in% class(x))) {
      x <- format(x)
      x <- substring(x, 2, nchar(x)-4)
      data[[i]] <- x
    }
  }

  w <- latex(data, file=file.path(TexDirName(), sprintf("%s.tex", panel)),
             title=panel, colheads=lab,
             longtable=longtable, size=size, caption=caption,
             landscape=landscape, rowname=NULL)
  invisible()
}
