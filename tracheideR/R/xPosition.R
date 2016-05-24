xPosition = function (a, scale = 1) {
  x <- a$cut.points[,1] - a$cut.points[1,1]
  pos <- x[a$what == 2]
  list("RingWidth" = x[length(x)] * scale,
       "x"         = (pos[-c(1,length(pos))] - a$LD / 2) * scale)
}
