#' Cross-validation for spatio-temporal imputation
#' 
#' @param data a data matrix with missing values
#' @param FUN the imputation function to be evaluated by cross-validation
#' @param date.info logical, if date information is provided in the data. 
#' @param cfold fold size on the columns
#' @param rfold fold size on the rows.
#' @param ... other arguments
#' @return the cross-validated RMSE
#' 
#' @export
#' @examples
#' data(hqmr.data)
#' # the real cross-validation will take some time to finish
#' # impCV(hqmr.data)
#' 
impCV <- function (data, FUN = Cut, date.info = TRUE, 
                   cfold = 10, rfold = 10, ...) {
  FUN <- match.fun(FUN)
  if(!date.info){
  data <- as.matrix(data)
  } else {
      date <- data[, "date"]
      data <- as.matrix(subset(data, select = -date))
  }
  
  nr <- nrow(data)
  nc <- ncol(data)
  r <- sample(1:nr, nr)
  s <- sample(1:nc, nc)
  cv.data <- data.frame(data[r, s])
  
  if(date.info){
  cv.data$date <- date[r]
  }
  
  by_i <- ceiling(nr/rfold)
  by_j <- ceiling(nc/cfold)
  i_seq <- seq(1, nr, by = by_i)
  j_seq <- seq(1, nc, by = by_j)
  cv.rmse <- matrix(0, length(i_seq), length(j_seq))
  
  for (j in j_seq) {
    for (i in i_seq) {
      m <- (i %/% by_i) + 1
      n <- (j %/% by_j) + 1
      if (n <= (cfold -1)){
        grid <- cv.data[i:(i + (by_i - 1)), j:(j + (by_j - 1))]
        grid.old <- grid # xtrue
        cv.na <- cv.data
        cv.na[i:(i + (by_i - 1)), j:(j + (by_j - 1))] <- NA # 
        cv.new <- FUN(cv.na, ...) # 
        grid.new <- cv.new[i:(i + (by_i - 1)), j:(j + (by_j - 1))] # ximp
        cv.rmse[m, n] <- Grmse(grid.new, grid.old)
      }
      else {
        grid <- cv.data[i:(i + (by_i - 1)), j:(j + nc - j_seq[cfold])]
        grid.old <- grid
        cv.na <- cv.data
        cv.na[i:(i + (by_i - 1)), j:(j + nc - j_seq[cfold])] <- NA
        cv.new <- FUN(cv.na, ...)
        grid.new <- cv.new[i:(i + (by_i - 1)), j:(j + nc - j_seq[cfold])]
        cv.rmse[m, n] <- Grmse(grid.new, grid.old)
      }
    }
  }
  return(rmse = as.vector(cv.rmse))
}
