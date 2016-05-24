jll.gNegBinomial <- function(y, y.hat, groups, ...) {
  y.hat.sum.g <- tapply(y.hat, groups, sum)
  y.sum.g <- tapply(y, groups, sum)
  both.sum.g <- y.hat.sum.g + y.sum.g
  lg.y.hat.sum.g <- tapply(lgamma(y.hat), groups, sum)
  lg.y.sum.g <- tapply(lgamma(y + 1), groups, sum)
  lg.both.sum.g <- tapply(lgamma(y + y.hat), groups, sum)
  ll <- sum(lgamma(y.hat.sum.g) + lgamma(y.sum.g + 1) - 
            lgamma(both.sum.g) + lg.both.sum.g -
            lg.y.sum.g - lg.y.hat.sum.g)

}
