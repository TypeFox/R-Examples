write.wcmd <- function(cmd, filename) {

  out <- list(start = "&INST")

  cmdnames <- unique(c("title", "data","item1", "ni", "name1",
    "namelen", "csv", "hlines"), names(cmd)[!names(cmd) %in%
    c("codes", "tfile", "arglist", "anchor", "labels", "extra")])
  for(i in cmdnames) {
      if(!is.null(cmd[[i]]))
      out[[i]] <- paste(i, cmd[[i]], sep = "=")
  }

  if(!is.null(cmd$codes))
    out$codes <- paste("codes=",
      paste(cmd$codes, collapse = ""))

  if(!is.null(cmd$tfile))
    out$tfile <- c("TFILE=*", unlist(cmd$tfile), "*")

  for(i in names(cmd$arglist))
    out[[i]] <- paste(i, cmd$arglist[[i]], sep = "=")

  if(!is.null(cmd$anchor))
    out$anchor <- c("IAFILE=*", paste(cmd$anchor[, 1],
      cmd$anchor[, 2]), "*")

  out$extra <- cmd$extra

  out$end <- "&END"

  if(!is.null(labels))
    out$labels <- c(cmd$labels, "END LABELS")

  writeLines(unlist(out), con = filename)
}
