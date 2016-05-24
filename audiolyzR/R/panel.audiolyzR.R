panel.audio <- function(x, y, ...) {
  audioScatter (x = "x", y = "y",
                data = data.frame (x = x, y = y),
                show.plots = FALSE)
}

panel.audiolyzR <- function(x, y, ...){
  panel.audio (x, y, ...)
  panel.hexbinplot (x, y, ...)
}