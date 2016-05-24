mts2ts <- 
function(x, seas = 1:frequency(x), na.rm = FALSE) {
  if (!is.mts(x)) 
    stop("x must be of class 'mts'")
  st <- start(x)[1]
  x1 <- window(x, start = st, end = c(end(x)[1], frequency(x)), 
               extend = TRUE)
  x1 <- aggregate(x1, 1, meanSub, sub = seas, na.rm = na.rm)
  x1 <- as.numeric(t(x1))
  ts(x1, start = st, frequency = ncol(x))
} 
