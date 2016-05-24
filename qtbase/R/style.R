
## customizing appearance through stylesheets

construct.stylesheet <- function(..., what = "*")
{
    elts <- list(...)
    nms <-
        if (is.null(names(elts))) rep("", length(elts))
        else names(elts)
    ans <- character(length(elts))
    for (i in seq_along(elts))
    {
        ans[i] <-
            if (nms[i] == "") sprintf(" %s { %s } ", what, elts[i])
            else sprintf(" %s { %s : %s } ", what, nms[i], elts[i])
    }
    ans
}

qsetStyleSheet <- function(..., what = "*", widget = NULL, append = TRUE)
{
    style <- construct.stylesheet(..., what = "*")
    if (append) style <- c(qstyleSheet(widget), style)
    style <- paste(style, collapse = "\n")
    if (is.null(widget)) {
      app <- Qt$QApplication$instance()
      app$styleSheet <- style
    } else {
      widget$styleSheet <- style
    }
}

qstyleSheet <- function(widget = NULL)
{
  if (is.null(widget)) {
    Qt$QApplication$instance()$styleSheet
  } else {
    widget$styleSheet
  }
}

