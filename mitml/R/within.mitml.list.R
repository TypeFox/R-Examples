within.mitml.list <- function(data, expr, ...){
# evaluate an expression for a list of data sets, then return altered data sets

  expr <- substitute(expr)
  parent <- parent.frame()

  out <- lapply(data, function(x){
    e <- evalq(environment(), x, parent)
    eval(expr, e)
    l <- as.list(e)
    l <- l[!sapply(l, is.null)]
    nD <- length(del <- setdiff(names(x), (nl <- names(l))))
    x[nl] <- l
    if(nD){
      x[del] <- if(nD==1){ NULL } else { vector("list", nD) }
    }
    x
  })
  class(out) <- c("mitml.list","list")
  out

}
