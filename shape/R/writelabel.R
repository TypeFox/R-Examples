
##==============================================================================
## writelabel   : adds a label next to a plot
##==============================================================================

writelabel <- function (text=NULL, nr=1, at=-0.1, line=1, cex=1.5, ...) {

  if (is.null(text))
    text <- LETTERS[nr]

  ## scale factors
  usr    <- par("usr")
  xmin   <- usr[1]
  xmax   <- usr[2]
  xrange <- xmax-xmin
  pos    <- xmin  + at * xrange

  mtext ( text=text, at=pos, line=line, cex=cex, ...)

}
