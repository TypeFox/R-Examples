organize.terms <-
function(obj) {
  if (obj$recursive) {
    children.copy <- obj$children
    obj$children <- list()
    rec <- unlist(lapply(children.copy, FUN = function(x) x$recursive))
    children.rec <- children.copy[rec]
    if (length(children.rec) > 0) for(i in 1:length(children.rec)) obj$children[[i]] <- organize.terms(children.rec[[i]])    
    children.nonrec <- children.copy[!rec]
    if (length(children.nonrec) > 0) {
      ord <- order(unlist(lapply(children.nonrec, FUN = function(x) length(x$cond))), decreasing = TRUE)
      obj$children <- c(children.nonrec[ord], obj$children)  
    }
    if (obj$fraction) obj$divisor <- organize.terms(obj$divisor)
  }
  return(obj)
}
