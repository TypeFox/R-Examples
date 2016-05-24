`seasonal.ave` <-
function(x, ave.FUN = mean, ...)
{
  result <- ts(start = c(0, 1), end = c(1, 0),
               frequency = frequency(x))
  result[cycle(x)] <- ave(x, cycle(x), FUN = ave.FUN, ...)
  return(result)
}
