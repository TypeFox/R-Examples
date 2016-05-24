is.exTabular <- function(x) {
  is(x, "exTabular")
}

#' @export
print.exTabular <- function (x, ...) {
  print(x$tables)
}

#' @export
summary.exTabular <- function (object, ...) {
  print(object$tags$title)
  print("---------------------------------")
  print(object$tables)
}

#Anonymous constructor
.exTabular <- function(tables, target = NULL, formats, tableSplit, tableType, title, tags, alias="table") {
  
  if (!is.null(target))
    newTags <- .metaTags(title = title, target = target, alias = alias)
  else
    newTags <- .metaTags(title = title, alias = alias)
  tags <- .updateTags(tags, newTags)
  
  tab <- list(
    "tables"     = tables,#a list of data.frames (one per output) with the data
    "formats"    = formats,#a list of data.frames (one per output, each one same dim of data) with the format and appearance of data
    "tableSplit" = tableSplit,
    "tableType"  = tableType,
    "tags"       = tags
  )
  
  class(tab) <- c("exTabular", "reportable")
  
  tab
}