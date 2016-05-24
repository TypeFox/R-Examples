mokkenTable <-
function (data) {

  a <- data.frame(summary(mokken::check.monotonicity(stats::na.omit(data))))[,c(1,9,10)]
  b <- data.frame(summary(mokken::check.iio(stats::na.omit(data)))$item.summary)[,c(9,10)]
  table <- cbind(a, b)
  
  names(table) <- c("Hi", "#mono-sig", "mono-crit", "#iio-sig", "iio-crit")
  table
}
