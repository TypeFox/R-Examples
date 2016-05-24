# based on ideas by Miquel-Angel Garcia-Lopez
gwr.tricube <- function (dist2, d) {
      d3 <- d^3
      dist3 <- dist2^1.5
      w <- ifelse(dist3 > d3, 0, (1-(dist3/d3))^3)
      w
}

