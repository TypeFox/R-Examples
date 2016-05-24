hist.lengths <-
function (x, ..., log = FALSE, zeros.rm = TRUE) {
  newlen <- x$length + 0
  if (zeros.rm & x$zeros) {
    newlen <- x$length + x$maxcens
  }
  args <- list(...)
  breaks <- "Sturges"
  if (!is.null(args$breaks)) breaks <- args$breaks
  include.lowest <- TRUE
  if (!is.null(args$include.lowest)) include.lowest <- args$include.lowest
  right <- TRUE
  if (!is.null(args$right)) right <- args$right
  res <- tapply(newlen, x$categories, function(Lengths) {
           if (log) {
             LogLengths <- log(Lengths)
             hist(LogLengths, breaks = breaks, include.lowest = include.lowest,
                  right = right, plot = FALSE)
           }
           else {
             hist(Lengths, breaks = breaks, include.lowest = include.lowest,
                  right = right, plot = FALSE)
           }
         })
  res$direction <- x$direction
  res$log <- log
  class(res) <- "hist.lengths"
  plot <- TRUE
  if (!is.null(args$plot)) plot <- args$plot
  if (plot) {
    plot(res, ...)
    return(invisible(res))
  }
  if (!plot) return(res)
}

