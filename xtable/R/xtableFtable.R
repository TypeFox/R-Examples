### ftable objects, requested by Charles Roosen
### Feature request #2248, 2/9/2012
xtableFtable <- function(x, caption = NULL, label = NULL, align = NULL,
                         digits = 0, display = NULL,
                         quote = FALSE,
                         method = c("non.compact", "row.compact",
                                    "col.compact", "compact"),
                         lsep = " $\\vert$ ", ...) {
  method <- match.arg(method)
  saveMethod <- method
  xDim <- dim(x)
  nRowVars <- length(attr(x, "row.vars"))
  nColVars <- length(attr(x, "col.vars"))
  if (nRowVars == 0){
    if (method =="col.compact"){
      method <- "non.compact"
    } else if (method == "compact"){
      method <- "row.compact"
    }
  }
  if (nColVars == 0){
    if (method =="row.compact"){
      method <- "non.compact"
    } else if (method == "compact"){
      method <- "col.compact"
    }
  }
  if (method == "non.compact"){
    nCharCols <- nRowVars + 2
    nCharRows <- nColVars + 1
  }
  if (method == "row.compact"){
    nCharCols <- nRowVars + 2
    nCharRows <- nColVars
  }
  if (method == "col.compact"){
    nCharCols <- nRowVars + 1
    nCharRows <- nColVars + 1
  }
  if (method == "compact"){
    nCharCols <- nRowVars + 1
    nCharRows <- nColVars
  }

  if(is.null(align)) {
    align <- c(rep("l", nCharCols - 1), "l |", rep("r", xDim[2]))
  }
  if(is.null(display)) {
    display <- c(rep("s", nCharCols), rep("d", xDim[2]))
  }

  attr(x, "ftableCaption") <- caption
  attr(x, "ftableLabel") <- label
  attr(x, "ftableAlign") <- align
  attr(x, "ftableDigits") <- digits
  attr(x, "quote") <- quote
  attr(x, "ftableDisplay") <- display
  attr(x, "method") <- method
  attr(x, "lsep") <- lsep
  attr(x, "nChars") <- c(nCharRows, nCharCols)
  class(x) <- c("xtableFtable", "ftable")
  return(x)
}

print.xtableFtable <- function(x,
  type = getOption("xtable.type", "latex"),
  file = getOption("xtable.file", ""),
  append = getOption("xtable.append", FALSE),
  floating = getOption("xtable.floating", TRUE),
  floating.environment = getOption("xtable.floating.environment", "table"),
  table.placement = getOption("xtable.table.placement", "ht"),
  caption.placement = getOption("xtable.caption.placement", "bottom"),
  caption.width = getOption("xtable.caption.width", NULL),
  latex.environments = getOption("xtable.latex.environments", c("center")),
  tabular.environment = getOption("xtable.tabular.environment", "tabular"),
  size = getOption("xtable.size", NULL),
  hline.after = getOption("xtable.hline.after", NULL),
  NA.string = getOption("xtable.NA.string", ""),
  only.contents = getOption("xtable.only.contents", FALSE),
  add.to.row = getOption("xtable.add.to.row", NULL),
  sanitize.text.function = getOption("xtable.sanitize.text.function", as.is),
  sanitize.rownames.function = getOption("xtable.sanitize.rownames.function",
                                         sanitize.text.function),
  sanitize.colnames.function = getOption("xtable.sanitize.colnames.function",
                                         sanitize.text.function),
  math.style.negative = getOption("xtable.math.style.negative", FALSE),
  math.style.exponents = getOption("xtable.math.style.exponents", FALSE),
  html.table.attributes = getOption("xtable.html.table.attributes", "border=1"),
  print.results = getOption("xtable.print.results", TRUE),
  format.args = getOption("xtable.format.args", NULL),
  rotate.rownames = getOption("xtable.rotate.rownames", FALSE),
  rotate.colnames = getOption("xtable.rotate.colnames", FALSE),
  booktabs = getOption("xtable.booktabs", FALSE),
  scalebox = getOption("xtable.scalebox", NULL),
  width = getOption("xtable.width", NULL),
  comment = getOption("xtable.comment", TRUE),
  timestamp = getOption("xtable.timestamp", date()),
  ...) {
  if (type == "latex"){
    ## extract the information in the attributes
    caption <- attr(x, "ftableCaption")
    label <- attr(x, "ftableLabel")
    align <- attr(x, "ftableAlign")
    digits <- attr(x, "ftableDigits")
    quote <- attr(x, "quote")
    digits <- attr(x, "ftabelDigits")
    method <- attr(x, "method")
    lsep <- attr(x, "lsep")
    nCharRows <- attr(x, "nChars")[1]
    nCharCols <- attr(x, "nChars")[2]
    nRowVars <- length(attr(x, "row.vars"))
    nColVars <- length(attr(x, "col.vars"))
    
    ## change class so format method will find format.ftable
    ## even though format.ftable is not exported from 'stats'
    class(x) <- "ftable"
    fmtFtbl <- format(x, quote = quote, digits = digits,
                      method = method, lsep = lsep)
    attr(fmtFtbl, "caption") <- caption
    attr(fmtFtbl, "label") <- label

    ## sanitization is possible for row names and/or column names
    ## row names
    if (is.null(sanitize.rownames.function)) {
      fmtFtbl[nCharRows, 1:nRowVars] <-
        sanitize(fmtFtbl[nCharRows, 1:nRowVars], type = type)
    } else {
      fmtFtbl[nCharRows, 1:nRowVars] <-
        sanitize.rownames.function(fmtFtbl[nCharRows, 1:nRowVars])
    }
    ## column names
    if (is.null(sanitize.colnames.function)) {
      fmtFtbl[1:nColVars, nCharCols - 1] <-
        sanitize(fmtFtbl[1:nColVars, nCharCols - 1],
                 type = type)
    } else {
      fmtFtbl[1:nColVars, nCharCols - 1] <-
        sanitize.colnames.function(fmtFtbl[1:nColVars, nCharCols - 1])
    }
    ## rotations are possible
    if (rotate.rownames){
      fmtFtbl[1:dim(fmtFtbl)[1], 1:(nCharCols - 1)] <-
        paste0("\\begin{sideways} ",
               fmtFtbl[1:dim(fmtFtbl)[1], 1:(nCharCols - 1)],
               "\\end{sideways}")
    }
    if (rotate.colnames){
      if (rotate.rownames){
        fmtFtbl[1:(nCharRows), (nCharCols):dim(fmtFtbl)[2]] <-
          paste0("\\begin{sideways} ",
                 fmtFtbl[1:(nCharRows), (nCharCols):dim(fmtFtbl)[2]],
                 "\\end{sideways}")
      } else {
        fmtFtbl[1:(nCharRows), 1:dim(fmtFtbl)[2]] <-
          paste0("\\begin{sideways} ",
                 fmtFtbl[1:(nCharRows), 1:dim(fmtFtbl)[2]],
                 "\\end{sideways}")
      }
    }


    ## booktabs is incompatible with vertical lines in tables
    if (booktabs) align <- gsub("|","", align, fixed = TRUE)
    attr(fmtFtbl, "align") <- align
    attr(fmtFtbl, "digits") <- digits
    attr(fmtFtbl, "quote") <- quote
    attr(fmtFtbl, "display") <- display

    ## labels should be left aligned
    for (i in 1:nCharRows){
      fmtFtbl[i, nCharCols:dim(fmtFtbl)[2]] <-
        paste0("\\multicolumn{1}{l}{ ",
               fmtFtbl[i, nCharCols:dim(fmtFtbl)[2]], "}")
    }


    if(is.null(hline.after)) {
      hline.after <- c(-1, nCharRows, dim(fmtFtbl)[1])
    }
    print.xtable(fmtFtbl, hline.after = hline.after,
                 include.rownames = FALSE, include.colnames = FALSE,
                 booktabs = booktabs,
                 sanitize.text.function = as.is)
  } else {
    stop("print.xtableFtable not yet implemented for this type")
  }
}

## format.xtableFtable <- function(x, quote = TRUE, digits = getOption("digits"),
##                                 method = c("non.compact", "row.compact",
##                                            "col.compact", "compact"),
##                                 lsep = " | ", ...){
##   class(x) <- "ftable"
  
##   format(x, quote = quote, digits = digits,
##          method = method, lsep = lsep, ...)
## }
