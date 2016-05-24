.parseRVignetteMetadata <- function(text, ...) {
  # Parse "\Vignette" directives into RSP metadata
  bfr <- unlist(strsplit(text, split="\n", fixed=TRUE), use.names=FALSE)

  pattern <- "[[:space:]]*%*[[:space:]]*\\\\Vignette(.*)\\{([^}]*)\\}"
  keep <- (regexpr(pattern, bfr) != -1L)
  bfr <- bfr[keep]

  # Nothing todo?
  if (length(bfr) == 0L) return(list())

  # Mapping from R vignette metadata to RSP metadata
  map <- c(
    # Official R vignette markup
    "IndexEntry"="title",
    "Keyword"="keyword",
    "Keywords"="keywords", ## Deprecated in R
    "Engine"="engine",
    # Custom
    "Subject"="subject",
    "Author"="author",
    "Date"="date",
    "Tangle"="tangle",
    "Compression"="compression"
  )

  metadata <- grep(pattern, bfr, value=TRUE)
  names <- gsub(pattern, "\\1", metadata)
  metadata <- gsub(pattern, "\\2", metadata)
  metadata <- trim(metadata)

  # Keep only known markup
  keep <- is.element(names, names(map))
  metadata <- metadata[keep]
  names <- names[keep]

  # Nothing todo?
  if (length(names) == 0L) return(list())

  # Rename
  names <- map[names]
  names(metadata) <- names
  metadata <- as.list(metadata)

  # Special: Merge all keyword meta data into one comma-separated entry
  idxs <- which(is.element(names(metadata), c("keyword", "keywords")))
  keywords <- unlist(metadata[idxs], use.names=FALSE)
  keywords <- paste(keywords, collapse=", ")
  metadata <- metadata[-idxs]
  metadata$keywords <- keywords

  metadata
} # .parseRVignetteMetadata()


##############################################################################
# HISTORY:
# 2015-05-11
# o Added \VignetteCompression{}.
# 2015-02-04
# o Now all \VignetteKeyword{} directives are parsed into meta data
#   field 'keywords'.
# 2013-09-18
# o Now .parseRVignetteMetadata() also records R vignette meta data
#   'engine' (from %\VignetteEngine{}) and RSP custom 'tangle' (from
#   %\VignetteTangle{}).
# o Added .parseRVignetteMetadata().  Was an internal function of the
#   preprocess() method for RspDocument.
##############################################################################
