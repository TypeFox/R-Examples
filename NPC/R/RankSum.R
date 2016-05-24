RankSum <- function (y, tr, tl, ...) {
  ## Rank-sum
  sum(rank(y, na.last=FALSE)[tr==tl])
}
