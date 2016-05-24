addpi <-
function(k, global=0.5, ...){
  x = seq(0, 1, 0.01)
  y = sapply(x, function(pr) piFunction(pr, kparam=k, qparam=global))
  lines(x, y, ...)
}
