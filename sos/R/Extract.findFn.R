"[.findFn" <- function(x, i, j,
    drop = if (missing(i)) TRUE else length(cols) == 1) {
##
## 1.  missing(j)
##
  class(x) <- 'data.frame'
  if(missing(j)){
    xi <- x[i, ]
    attr(xi, 'PackageSummary') <- PackageSummary(xi)
    class(xi) <- c("findFn", 'data.frame')
    return(xi)
  }
##
## 2.  select columns, return data.frame
##
  cols <- j
  x[i, j, drop=drop]
}
