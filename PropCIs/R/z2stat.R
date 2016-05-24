z2stat <-
function (p1x,nx,p1y,ny,dif){

      diff = p1x-p1y-dif
      if ( abs(diff) == 0 ) {
        fmdiff = 0}
      else{
        t = ny/nx
        a = 1+t
        b = -(1+ t + p1x + t*p1y + dif*(t+2))
        c = dif*dif + dif*(2*p1x + t +1) + p1x + t*p1y
        d = -p1x*dif*(1+dif)
        v = (b/a/3)^3 - b*c/(6*a*a) + d/a/2
        s = sqrt( (b/a/3)^2 - c/a/3)
        if(v>0){u=s}
         else{u=-s}
        w = (3.141592654+acos(v/u^3))/3
        p1d = 2*u*cos(w) - b/a/3
        p2d = p1d - dif
        nxy = nx + ny
        var = (p1d*(1-p1d)/nx + p2d*(1-p2d)/ny) * nxy / (nxy - 1) ## added: * nxy / (nxy - 1)
        fmdiff = diff^2/var
      }
  return(fmdiff)
}

