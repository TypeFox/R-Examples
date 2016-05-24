`rotation` <-
function(coords,angle)
{
radian <- (angle * pi)/180

x <- c(cos(radian),- sin(radian))
y <- c(sin(radian),cos(radian))
nlecoord=coords %*% cbind(x,y)

return(nlecoord)
  }

