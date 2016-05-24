e.bars <-
function(graph, m, ebl, sides=2, length=0) {
  if(sides==1) {
    l <- ifelse(m <= 0, -ebl, ebl)
    arrows(graph, m, graph, m+l, angle=90, code=3, length=length)
  }
  if(sides==2) {
    arrows(graph, m + ebl, graph, m - ebl, angle=90, code=3, length=length)
  }
}
