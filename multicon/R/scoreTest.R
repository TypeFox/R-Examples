scoreTest <- function(items, keys, Zitems = FALSE, maxScore = NULL, minScore = NULL, rel = FALSE, nomiss = .8, tr = 0, item.names = NULL, check.keys = TRUE) {
  if(is.null(maxScore)) {
    maxScore <- max(items, na.rm=TRUE)
  }
  if(is.null(minScore)) {
    minScore <- min(items, na.rm=TRUE)
  }
  if(Zitems) {
    items <- data.frame(scale(items))
    L <- lapply(keys, function(x) matrix(unlist(ifelse(x < 0, -1 * items[,abs(x)], items[,abs(x)])), ncol=length(x), byrow=FALSE))
  }
  else {
    L <- lapply(keys, function(x) matrix(unlist(ifelse(x < 0, (maxScore+minScore) - items[,abs(x)], items[,abs(x)])), ncol=length(x), byrow=FALSE))
  }
  scores.list <- lapply(L, function(x) apply(x, 1, meanif, nomiss=nomiss, tr=tr))
  if(is.null(item.names)) {
    item.names <- names(keys)
  }
  if(rel) {
    rel.list <- lapply(L, function(x) alpha(data.frame(x), check.keys=check.keys)$total)
    rels <- data.frame(matrix(unlist(rel.list), nrow=length(rel.list), byrow=TRUE))
    colnames(rels) <- c('raw_alpha', 'std.alpha', 'G6(smc)', 'average_r', 'S/N', 'ase', 'mean', 'sd')
    rownames(rels) <- item.names
    scores <- data.frame(matrix(unlist(scores.list), ncol=length(scores.list), byrow=FALSE))
    colnames(scores) <- item.names
    out <- list("rel"=rels, "scores"=scores)
  }
  else {
    out <- data.frame(matrix(unlist(scores.list), ncol=length(scores.list), byrow=FALSE))
    colnames(out) <- item.names
  }
  out
}