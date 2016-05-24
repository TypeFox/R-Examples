"chessboard" <-
function (color="black", seed=NULL)
{
  if(length(seed)) set.seed(seed)
  canvas()
  eighths <- seq(0, 1, length=9)
  boundary <- c(0, rep(rep(0:1, 8), each=2), 1)
  doub <- rep(eighths, each=2)
  doubx <- doub + c(0,0, runif(14, -.12, .12), 0, 0)
  douby <- doub + c(0,0, runif(14, -.12, .12), 0, 0)
  x <- c(doubx, 1 - boundary)
  y <- c(boundary, rev(douby))
  polygon(x[-18:-35], y[-18:-35], col=color)
}

