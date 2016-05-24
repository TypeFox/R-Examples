StereoCirc <-
function(n.seg = 360){
  cir.x <- cos(seq(from = 0, to  = 2 * pi, length = n.seg))
  cir.y <- sin(seq(from = 0, to  = 2 * pi, length = n.seg))
  lines(cir.x, cir.y)
  lines(c(-0.025, 0.025), c(0, 0), lwd = .5)
  lines(c(0, 0), c(-0.025, 0.025), lwd = .5)
}
