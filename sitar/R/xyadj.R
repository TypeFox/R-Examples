xyadj <- function(x, y=NULL, id, object, abc=ranef(object)[id, ], tomean=TRUE) {
#	returns x and y adjusted for random effects a, b and c
  abc[, letters[1:3][!letters[1:3] %in% names(ranef(object))]] <- 0 # omit not in model
  xoffset <- object$xoffset
  if (!is.na(b0 <- fixef(object)['b'])) xoffset <- xoffset + b0
  if (tomean) {
    x.adj <- (x - xoffset - abc$b) * exp(abc$c) + xoffset
    y.adj <- y - abc$a
  } else {
    x.adj <- (x - xoffset) / exp(abc$c) + xoffset + abc$b
    y.adj <- y + abc$a
  }
  return(list(x=x.adj, y=y.adj))
}
