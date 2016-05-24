## author: Olivia Lau
## date  : Dec, 2007
convert <- function(vars, data, order, ...) {
  vid <- unique(c(as.character(data[[match(vars, names(data))]])))
  if (!identical(sort(na.omit(order)), sort(na.omit(vid)))) {
    stop("all values must be identified in order")
  }
  if (NA %in% order) keepna <- NULL
  else keepna <- NA
  for (i in vars) {
    tmp <- data[[i]]
    tmp <- unclass(factor(tmp, levels = order, exclude = keepna, ...))
    data[[i]] <- tmp
   attr(data[[i]], "levels") <- NULL
  }
  data
}
