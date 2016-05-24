download.lib <- function(pkgs, destdir = getwd(), mode = 'wb',
  pdf.url = "http://cran.r-project.org/web/packages/", f.zip = TRUE,
  f.pdf = TRUE) 
{
  pa <- data.frame(available.packages(type="source"), stringsAsFactors=F)
  if (!all(pkgs %in% pa[, "Package"])) {
    stop("Some package names are incorrect.\n")}  
  pb <- pa[pkgs, , drop = FALSE]

  for( i in seq_along(pkgs)) {
    k <- pkgs[i]
    n <- paste(k, "_", pb[i, "Version"], ".tar.gz", sep = "")
 
    download.file(url = paste(pb[i, "Repository"], n, sep="/"),
      destfile = paste(destdir, n, sep="/"))
    if (f.zip) {download.packages(pkgs = pkgs[i], destdir = destdir)}
    if (f.pdf) {
      download.file(url = paste(pdf.url, k, "/", k, ".pdf", sep=""),
        destfile = paste(destdir, "/", k, ".pdf", sep=""), mode = mode)
    }
  }
  invisible(pkgs)
} 