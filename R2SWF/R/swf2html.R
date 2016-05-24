#' Embed the SWF file into an HTML page
#'
#' This function will generate an HTML file to display the Flash animation.
#' @param swf.file the path of the SWF file
#' @param output the output path of the HTML file; by default \file{foo.swf}
#'   produces \code{foo.html} if not specified (set \code{FALSE} so that no file
#'   will be written)
#' @param width width of the Flash
#' @param height height of the Flash
#' @param fragment whether to produce an HTML fragment only
#' @return The HTML code as a character string.
#' @export
#' @author Yihui Xie <\url{http://yihui.name}>
#' @examples
#' output = dev2swf({
#'   for (i in 1:10) plot(runif(20), ylim = c(0, 1))
#' }, output = 'test.swf')
#' swf2html(output)
swf2html = function(swf.file, output, width = 480, height = 480, fragment = FALSE) {
  if (!file.exists(swf.file)) stop("swf file does not exist")
  if (missing(output)) output = basename(sub('\\.swf$', '.html', swf.file))
  size = paste(c(sprintf('width="%s"', width), sprintf('height="%s"', height)), collapse = ' ')
  html = sprintf('<embed %s name="plugin" src="%s" type="application/x-shockwave-flash">',
                 size, swf.file)
  if (!fragment) html = paste('<html>
<head>
  <title>Flash animations made by the R2SWF package</title>
</head>
<body>
<div align="center">
', html, '
</div>
</body>
</html>
')
  if (!identical(output, FALSE)) cat(html, file = output)
  if (is.character(output) && file.exists(output)) {
      if (dirname(output) != '.') {
          file.rename(output, basename(output))
          output = basename(output)
      }
    message('output at ', normalizePath(output))
    if (interactive()) try(browseURL(output), silent = TRUE)
  }
  invisible(html)
}
