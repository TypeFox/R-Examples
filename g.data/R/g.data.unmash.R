g.data.unmash <- function(fn) gsub("@", "", sub("\\.RData$", "", basename(fn)))
