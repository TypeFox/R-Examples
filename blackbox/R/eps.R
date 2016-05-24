## a user-friendly version of postscript(), see http://r.789695.n4.nabble.com/Defaults-for-postscript-td918013.html
eps <- function(file = "Rplots.eps", height = NULL,
               width =  NULL, ...) {
  if (is.null(height)) height <- blackbox.getOption("epsArgs")$height
  if (is.null(height)) height <- 7.5 ## 6 not enough for rawProfiles at least ('figure margins too large')
  if (is.null(width)) width <- blackbox.getOption("epsArgs")$width
  if (is.null(width)) width <- 7.5
  postscript(file = file, onefile = TRUE, horizontal = FALSE,
             height = height, width = width, paper = "special", ...)
}
