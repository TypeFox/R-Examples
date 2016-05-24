.html.end <- function(file) {
  out <- '</body></html>'
  cat(out, file=file, append=TRUE)
}
