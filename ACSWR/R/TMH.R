TMH <-
function(x) {
  qh <- fivenum(x,c(0.25,0.5,0.75))
  return((qh[2]+(qh[1]+qh[3])/2)/2)
}
