skewdivTwo <-
function(xb1, xb2) {

  skew1 <- apply(xb1, 1, function(y) mean(y) - median(y))
  skew2 <- apply(xb2, 1, function(y) mean(y) - median(y))

  CDF1 <- ecdf(skew1)
  CDF2 <- ecdf(skew2)

  skew1sort <- sort(skew1)
  skew2sort <- sort(skew2)
  skewrange <- range(c(skew1, skew2))

  abs(sum(diff(c(skewrange[1], skew1sort, skewrange[2]))*CDF1(c(skewrange[1], skew1sort))) - 
    sum(diff(c(skewrange[1], skew2sort, skewrange[2]))*CDF2(c(skewrange[1], skew2sort))))

}
