harmonicspot <- function (z,x0,y0,cx,cy,b.x,b.y,retention,var.x,var.y) { 
  xt<-cx +b.x*cos(z)
  yt<-cy +b.y*cos(z)+retention*sin(z)
  dist<-sqrt(((xt-x0)^2)/var.x+((yt-y0)^2)/var.y)
  return("dist"=dist)}
