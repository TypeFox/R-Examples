read.wdat <- function(cmd, datfile = cmd$data, na = " ",
  ilabels = NULL) {

  dat <- readLines(datfile)
  ncols <- max(nchar(dat))
  nvars <- cmd$ni + 1

  if(!all(nchar(dat) == ncols))
    warning("Lines", paste(which(nchar(dat) != ncols),
      collapse = ", "), "are missing information")

  dat <- format(dat)
  out <- data.frame(matrix(ncol = nvars, nrow = length(dat)))
  out[, 1] <-
    trim(substr(dat, cmd$name1, cmd$name1 + cmd$namelen - 1))
  for(i in 1:cmd$ni) {
    out[, i + 1] <- substr(dat, cmd$item1 + i - 1,
      cmd$item1 + i - 1)
    out[out[, i + 1] == na, i + 1] <- NA
    out[, i + 1] <- as.numeric(out[, i + 1])
  }

  if(!is.null(ilabels))
    cmd$labels <- ilabels
  if(!is.null(cmd$labels))
    names(out) <- c("name", cmd$labels)
  else
    names(out) <- c("name", paste("i", 1:cmd$ni, sep = ""))

  return(out)
}
