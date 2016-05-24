imagesc <-
  function(G, col=grey((1:99)/100), ... )
  {
#########  plot an image after flipping and transposing
###   to match the matlab code
    d = dim(G)
    b = G[d[1]:1,]

    x = attr(G, "x")
    y = attr(G, "y")

    if(!is.null(x) & !is.null(y))
      {
        image(x=x, y=y, z=t(b), col=col, ...)
      }
    else
      {
        image(t(b), col=col, ...)
      }

    
    
    
  }
contoursc <- function (G, ...) 
{
    d = dim(G)
    b = G[d[1]:1, ]
    x = attr(G, "x")
    y = attr(G, "y")
    if (!is.null(x) & !is.null(y)) {
        contour(x = x, y = y, z = t(b),  ...)
    }
    else {
        contour(t(b),  ...)
    }
}
