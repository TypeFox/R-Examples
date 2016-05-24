ellipsespot <- function (z,x0,y0,cx,cy,semi.major,semi.minor,rote.rad) { 
  xt<-cx +semi.major*cos(rote.rad)*cos(z)-semi.minor*sin(rote.rad)*sin(z)
  yt<-cy +semi.major*sin(rote.rad)*cos(z)+semi.minor*cos(rote.rad)*sin(z)
  dist<-(xt-x0)^2+(yt-y0)^2
  return("dist"=dist)}
