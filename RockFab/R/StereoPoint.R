StereoPoint <-
function(my.az = 90, my.inc = 45, my.color = 'black', my.pch = 19, my.size = .25, my.label){
  my.az <-  my.az * (pi / 180) + pi
  my.inc <- my.inc * (pi / 180) - pi / 2
  my.tq <- sqrt(2) * sin(my.inc / 2)
  my.x <- my.tq * sin(my.az)
  my.y <- my.tq * cos(my.az)
  
  if(!missing(my.label)){
    i <- 0
    par(ps = 8)
    while(i < length(my.label)){
      i <- 1 + i
      text(my.x[i] + .025, my.y[i] + .025, my.label[i], cex = .9, adj = c(0, 0))
    }
  }
  points(my.x, my.y, pch = my.pch, col = my.color, cex = my.size)
}
