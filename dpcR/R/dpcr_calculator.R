dpcr_calculator <- function(k, n, average = FALSE, log = FALSE) {
  if (n < k)  
    stop("'n' must be larger or equal to 'k'.")
  if (k == 0)  
    stop("'k' must be larger then 0. No template molecules in the sample")
  if (!(average %in% c(TRUE, FALSE)))
    stop("'average' must be TRUE or FALSE.")
  
  sig5 <- round(sqrt(k*(1 - k/n))*5, 0)
  start_r <- ifelse(k - sig5 < 0, 0, k - sig5)
  stop_r <- ifelse(k + sig5 >= n, n - 1, k + sig5)
  range <- start_r:stop_r
  res <- dbinom(range, size = n, prob = k/n, log = log)
  if (average)
    range <- -log(1 - range/n)
  matrix(c(range, res), ncol = 2)
}