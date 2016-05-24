ct2df <- function(
  x,
  bid1  = "bid1",  
  bid2h = "bidh",  
  bid2l = "bidl",  
  yy    = "yy",
  yn    = "yn",
  ny    = "ny",
  nn    = "nn",
  y     = "y", 
  n     = "n", 
  type  = "double" 
)
{
  # single-bounded
  if (type == "single") {

    if (ncol(x) != 3) {
      stop("number of columns of x must be 3 for single-bounded")
    }

    reshape.x <- x[c(bid1, y, n)]

    bid1 <- c(rep(reshape.x[, 1], reshape.x[, 2]),
              rep(reshape.x[, 1], reshape.x[, 3])) 
    R1   <- c(rep(1L, sum(reshape.x[, 2])),
              rep(0L, sum(reshape.x[, 3])))

    cv.data <- data.frame(R1 = R1, bid1 = bid1)

  } else { # double-bounded

    if (ncol(x) != 7) {
      stop("number of columns of x must be 7 for double-bounded")
    }

    reshape.x <- x[c(bid1, bid2h, bid2l, yy, yn, ny, nn)]
    colnames(reshape.x) <- c("B1", "B2H", "B2L", "yy", "yn", "ny", "nn")

    bid.table      <- reshape.x[, c(1, 2, 3)]
    response.table <- reshape.x[, c(4, 5, 6, 7)]

    B1 <- rep(bid.table[, 1], rowSums(response.table))
    R  <- rep(names(response.table), response.table[1, ])
    for (i in 2:nrow(x)) {
      R <- c(R, rep(names(response.table), response.table[i, ]))
    }

    data <- data.frame(
      B1 = B1,
      R  = R,
      R1 = c(R == "yy") + c(R == "yn"),
      R2 = c(R == "yy") + c(R == "ny"))

    cv.data <- merge(bid.table, data,  by = "B1")
    cv.data$bid1 <- cv.data$B1
    cv.data$bid2 <- cv.data$B2H * cv.data$R1 + cv.data$B2L * (cv.data$R1 == 0)
  }

  return(cv.data)
}
