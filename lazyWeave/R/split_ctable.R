#' @rdname WritePrintCtable
#' @export split_ctable
#' 

split_ctable <- function(x, max.rows=35, keepVarTogether=TRUE, ...){
  lazy.args <- list(...)
  if ("counter" %in% names(lazy.args)) counter <- lazy.args$counter else counter <- "table"
  rowList <- list()
  rows <- 1:nrow(x)
  if (length(rows) < max.rows) rowList <- list(rows)
  else if (!keepVarTogether){
    ntab <- ceiling(nrow(x) / max.rows)
    index <- rep(1:ntab, rep(max.rows, ntab))[1:nrow(x)]
    rowList <- lapply(unique(index), function(x) rows[index==x])
  }
  else {
    while(length(rows) > max.rows){
      lastRow <- if (!is.na(x$name[rows[max.rows] + 1])) max.rows else max(which(!is.na(x$name[rows[1]:rows[max.rows]]))) - 1
      rowList <- c(rowList, list(rows[1]:rows[lastRow]))
      rows <- rows[-c(1:lastRow)]
      if (length(rows) < max.rows) rowList <- c(rowList, list(rows[1]:utils::tail(rows, 1)))
    }
  }
  Splits <- lapply(rowList, function(t) x[t, ])
  Tabs <- sapply(Splits, write.ctable, ...)
  Tabs[-1] <- paste("\\addtocounter{", counter, "}{-1}", Tabs[-1], sep="")
  paste(Tabs, collapse="\\clearpage")
}
