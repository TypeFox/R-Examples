write.wdat <- function(x, cmd, datfile = cmd$data, na = " ") {

  xnames <- colnames(x)

  if(is.null(xnames)) {
    i <- x
    p <- NULL
    warning("All columns in 'x' are unnamed and are assumed to ",
      "be item responses")
  }
  else {
    if("labels" %in% names(cmd))
      i <- x[, cmd$labels]
    else
      i <- x[, !xnames %in% "name"]
    if("name" %in% xnames)
      p <- x[, "name"]
    else
      p <- NULL
  }

  ncols <- max(cmd$item1 + cmd$ni, cmd$name1 + cmd$namelen) -
    cmd$namelen
  if(cmd$item1 > cmd$name1)
    iindex <- cmd$item1:(cmd$ni + cmd$item1 - 1) - cmd$namelen + 1
  else
    iindex <- cmd$item1:(cmd$ni + cmd$item1 - 1)

  out <- matrix(rep(" ", nrow(x)*ncols), ncol = ncols)
  if(is.null(p))
    iindex <- iindex - 1
  else
    out[, cmd$name1] <- format(as.character(p),
      width = cmd$namelen)
  out[, iindex] <- unlist(i)

  write.table(out, datfile, quote = FALSE, col.names = FALSE,
    row.names = FALSE, sep = "", na = na)
}
