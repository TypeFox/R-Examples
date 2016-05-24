##' @export
score.biprobit <- function(x,indiv=FALSE,...) {
  if (indiv) { s <- x$score; attributes(s)$logLik <- NULL; return(s) }
  colSums(x$score)
}
