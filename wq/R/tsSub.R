tsSub <- 
function(x1, seas = 1:frequency(x1)) {

  if (!is(x1, "ts"))
    stop("x1 must be of class 'ts'")
  stx <- start(x1)
  frx <- frequency(x1)
  if (!is(x1, "mts"))
    dim(x1) <- c(length(x1), 1)
  x2 <- window(x1, start = stx[1], end = c(end(x1)[1], frx), extend =
                 TRUE)
  x3 <- x2[cycle(x2) %in% seas, ]
  ts(x3, start = stx[1], frequency = length(seas))
}
