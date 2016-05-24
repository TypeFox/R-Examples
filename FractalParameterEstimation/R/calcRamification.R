calcRamification <-
function(figure) {
  count = 1
  while(figure > 3) {
    figure = figure/3
    count = increment(count)
  } # end while
  return(count)
}
