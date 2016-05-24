TestMatchDataFrameRowsOnly <- function() {
  x1 <- data.frame(larry = rnorm(10), moe = 1:10, curly = rpois(10, 2))
  x2 <- x1[c(1:5, 10:6), ]
  m <- MatchDataFrame(x1, x2)
  checkEqualsNumeric(m$column.permutation, 1:3)
  checkTrue(all(x2[m$row.permutation,
                   m$column.permutation] == x1))
  checkTrue(all(x2[m$row.permutation, ] == x1))
}

TestMatchDataFrameColumnsOnly <- function() {
  x1 <- data.frame(larry = rnorm(10), moe = 1:10, curly = rpois(10, 2))
  x2 <- x1[, c(3, 1, 2)]
  m <- MatchDataFrame(x1, x2)
  checkEqualsNumeric(m$row.permutation, 1:10)
  checkTrue(all(x2[m$row.permutation,
                   m$column.permutation] == x1))
  checkTrue(all(x2[, m$column.permutation] == x1))
}

TestMatchDataFrameRowsAndColumns <- function() {
  x1 <- data.frame(larry = rnorm(10), moe = 1:10, curly = rpois(10, 2))
  x2 <- x1[c(1:5, 10:6), c(3, 1, 2)]
  m <- MatchDataFrame(x1, x2)
  checkTrue(all(x2[m$row.permutation,
                   m$column.permutation] == x1))
}

TestNoMatch <- function() {
  checkStop(
      MatchDataFrame(data.frame(x = 1:3, y = 5:7),
                     data.frame(y = 6:8, x = 1:3)),
      expected = "Not all rows could be matched.",
      silent = TRUE)

  checkStop(
      MatchDataFrame(data.frame(x = 1:3, y = 5:7),
                     data.frame(z = 6:8, x = 1:3)),
      expected = "Not all columns could be matched.",
      silent = TRUE)

  checkStop(
      MatchDataFrame(data.frame(x = 1:3, y = 5:7),
                     data.frame(y = 6:9, x = 1:4)),
      expected = "Data are of different sizes",
      silent = TRUE)
}
