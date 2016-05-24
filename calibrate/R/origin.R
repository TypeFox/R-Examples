"origin" <-
function (m = c(0, 0), ...) 
{
   vx <- c(par("usr")[1],par("usr")[2])
   vy <- c(m[2],m[2])
   lines(vx,vy,...)
   vy <- c(par("usr")[3],par("usr")[4])
   vx <- c(m[1],m[1])
   lines(vx,vy,...)   
}
