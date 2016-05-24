#' @name ggplotGrid
#' @author Sven E. Templer
#' @title Arrange a List of ggplots
#' @description 
#' Arrange a list of ggplots with \link[gridExtra]{grid.arrange} and
#' output on local graphic device or as pdf/png when a path is supplied.
#' \code{ggplotGridA4} writes the plots to a DIN A4 (8 x 11 inches) pdf file
#' directly.
#' @param l List with ggplot objects.
#' @param path Plot to file of type pdf or png. Determine type by
#' path ending (.pdf or .png).\cr
#' Optional in \code{ggplotlist}: A character string that gives the path to export
#' the plot to a file, ending with 'pdf' or 'png' (case insensitive). If
#' missing, then the grid is returned to the current graphic device.
#' @param ncol Number of columns.
#' @param nrow Number of rows per page, only for pdfs.
#' @param width For pdfs/pngs the width in inches, else ignored.
#' @param height For pdfs/pngs the height in inches, else ignored.
#' @param res Resolution in dpi for pngs.
#' @param wide Wide format pdf pages (11x8 inches).
#' @param pdf.cairo Use \link{cairo_pdf} (or cairo_ps, svg) instead of \link{pdf}
#' @param onefile Create one file, see \link{cairo_pdf}.
#' @param ... Forwarded to cairo_pdf
#' @param x A list containing at least one ggplot object of class \code{gg}.
#' @examples
#' #
#' 
#' \dontrun{
#' library(ggplot2)
#' d <- data.frame(a=1:5,b=1:5)
#' x <- list(
#'   ggplot(d, aes(x=a,y=b,col=b)) + geom_line(),
#'   ggplot(d, aes(x=a,y=b,shape=factor(b))) + geom_point())
#' ggplotlist(x, 2)}
#' 
#' #

#' @rdname ggplotGrid
#' @export
ggplotGrid <- function (l, path, ncol = 1, nrow = 1,
                        width = 8, height = 11, res = 300,
                        pdf.cairo = TRUE, onefile = TRUE, ...) {
  
  # test the classes
  lclass <- sapply(l, function(i) class(i)[1] == "gg")
  if (!all(lclass))
    stop("Provide list with only ggplots!")
  
  # presets
  type <- 'int'
  n <- nrow * ncol # per pdf page
  ggempty <- list(
    ggplot(data.frame()) +
      geom_point() + theme_bw() +
      theme(panel.border = element_rect(color = "white")) +
      scale_x_discrete(breaks=NULL) +
      scale_y_discrete(breaks=NULL))
  
  # get device
  if (!missing(path)) {
    type <- strpart(path, "\\.", 123, roll=T)
    #type <- tolower(substr(path, nchar(path)-2, nchar(path)))
    if (type == 'pdf' && pdf.cairo)
      cairo_pdf(path, width = width, height = height, onefile = onefile, ...)
    else if (type == "ps")
      cairo_ps(path, width = width, height = height, onefile = onefile, ...)
    else if (type == "svg")
      svg(path, width = width, height = height, onefile = onefile, ...)
    else if (type == 'pdf')
      pdf(path, width = width, height = height)
    else if (type == 'png')
      png(path, width = width, height = height, units = 'in', res = res)
    else if (type == "eps") {
      setEPS()
      postscript(file = path, width = width, height = height, ...)
    } else
      stop(paste0('.', type, ' is an invalid ending, use: .pdf .png .eps .svg or .ps.'))
  }
  
  if (type == 'pdf') {
    avail <- TRUE
    while (avail) {
      ln <- length(l)
      # fill with empty to match grid
      if (ln < n) {
        empty <- (ln+1):n
        l[empty] <- ggempty
      }
      do.call(grid.arrange, c(l[1:n], ncol=ncol))
      l <- l[-c(1:n)]
      if (length(l) < 1) avail <- FALSE
    }
  } else {
    do.call(grid.arrange, c(l, ncol=ncol))
  }
  
  if (!missing(path)) dev.off()
  
  invisible(NULL)
  
}

#' @rdname ggplotGrid
#' @export
ggplotGridA4 <- function (l, path, ncol = 2, nrow = 1, wide = TRUE) {
  
  # test the classes
  lclass <- sapply(l, function(i) class(i)[1] == "gg")
  if (!all(lclass))
    stop("Provide list with only ggplots!")
  
  # set parameter
  n <- nrow * ncol
  width <- 11
  height <- 8
  if (!wide) {
    width <- 8
    height <- 11
  }
  
  # pdf
  pdf(path, width = width, height = height)
  avail <- TRUE
  while (avail) {
    ln <- length(l)
    # fill with empty to match grid
    if (ln < n) {
      empty <- (ln+1):n
      l[empty] <- list(ggplot(data.frame()) +
        geom_point() + theme_bw() +
        theme(panel.border = element_rect(color = "white")) +
        scale_x_discrete(breaks=NULL) +
        scale_y_discrete(breaks=NULL))
    }
    do.call(grid.arrange, c(l[1:n], ncol=ncol))
    l <- l[-c(1:n)]
    if (length(l) < 1) avail <- FALSE
  }
  dev.off()
  invisible(NULL)
  
}

#' @rdname ggplotGrid
#' @export ggplotlist
ggplotlist <- function (x, ncol = 1, path, width = 11, height = 8) {
  
  warning("ggplotlist is deprecated. Use 'ggplotGrid' or 'ggplotGridA4'. See the help page.")
  # check the classes
  keep <- sapply(x, function(i) inherits(i, "gg"))
  x <- x[keep]
  if (!any(keep))
    stop("Provide list with at least one ggplot.")
  else if (!all(keep))
    warning("Non-ggplot objects found and dropped.")
  
  # return
  if (!missing(path)) {
    type <- tolower(substr(path, nchar(path)-2, nchar(path)))
    if (type == 'pdf')
      pdf(path, width = width, height = height)
    else if (type == 'png')
      png(path, width = width, height = height, units = 'in', res = 300)
    else
      stop(paste0(type, ' is an invalid ending (use: pdf, png).'))
  }
  do.call(grid.arrange, c(x, ncol=ncol))
  if (!missing(path)) dev.off()
  invisible(NULL)
  
}
