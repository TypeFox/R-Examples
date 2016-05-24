search_files <- function(file, dir = ".", recursive = TRUE)
{
  if(is.null(dir)) dir <- "."
  lf <- list.files(dir, full.names = TRUE, recursive = recursive)
  bn <- basename(lf)
  gr <- lapply(file, function(f) which(bn == f))
  if(any((gl <- sapply(gr, length)) > 1L)) {
    warning(paste("The following files are found more than once, first match used:",
      paste(file[gl > 1L], collapse = ", ")))
  }
  ix <- rep.int(seq_along(file), gl)
  gr <- unlist(gr)
  fp <- file.path(dirname(lf[gr]), file[ix])
  fe <- file.exists(fp)
  fp[match(seq_along(file), ix[fe])]
}
