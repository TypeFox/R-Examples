"resid.squares" <-
function(x, y, y.hat, resid.plot="square", ...) {
  y.resid <- y - y.hat
  y.start <- ifelse(y.resid>0, y.hat, y)
  rect.height <- abs(y.resid)
  if (resid.plot=="square")
    rect.width <- rect.height *
      (par("pin")[2]/par("pin")[1]) /
        ((par("usr")[4]-par("usr")[3])/
         (par("usr")[2]-par("usr")[1]))
  else
    rect.width <- 0
  symbols(x+rect.width/2, y.start+rect.height/2,
          rectangles=cbind(rect.width,rect.height),
          inches=FALSE, add=TRUE, ...)
}
