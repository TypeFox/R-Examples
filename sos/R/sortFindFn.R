sortFindFn <- function(x, sortby=NULL) {
##
## 1.  PackageSummary
##
  x$Score <- as.numeric(as.character(x$Score))
  pkgSum <- PackageSummary(x, sortby)
##
## 2.  Sort order
##
  s0 <-  c('Count', 'MaxScore', 'TotalScore', 'Package',
           'Score', 'Function', 'Date', 'Description', 'Link')
  s0. <- tolower(s0)
  {
    if(is.null(sortby)) sortby <-  s0
    else {
      s1 <- match.arg(tolower(sortby), s0., TRUE)
      s1. <- c(s1, s0.[!(s0. %in% s1)])
      names(s0) <- s0.
      sortby <- s0[s1.]
    }
  }
##
## 3.  Merge(packageSum, x)
##
  packageSum <- pkgSum
  rownames(pkgSum) <- as.character(pkgSum$Package)
  pkgSum$Package <- NULL
  pkgS2 <- pkgSum[as.character(x$Package), , drop=FALSE]
  rownames(pkgS2) <- NULL
  xInS2 <- (names(x) %in% names(pkgS2))
  Ans <- cbind(as.data.frame(pkgS2), x[, !xInS2])
##
## 4.  Sort Ans by 'sort.'
##
  Ans.num <- Ans[, c('Count', 'MaxScore', 'TotalScore', 'Score')]
  ans.num <- cbind(as.matrix(Ans.num), Date=as.numeric(Ans$Date) )
  Ans.ch <- Ans[, c('Package','Function', 'Description', 'Link')]
  ans.ch <- as.data.frame(as.matrix(Ans.ch))
  ansKey <- cbind(as.data.frame(-ans.num), ans.ch)
#
  oSch <- do.call('order', ansKey[sortby])
  AnSort <- Ans[oSch,]
##
## 5.  attributes
##
  rownames(AnSort) <- NULL
#
#  attr(AnSort, "hits") <- hits
  attr(AnSort, 'PackageSummary') <- packageSum
#  attr(AnSort, 'string') <- string
#  attr(AnSort, "call") <- match.call()
  class(AnSort) <- c("findFn", "data.frame")
  AnSort
}
