#source('moonpos.R')
#source('sunpos.R')

mphase = function(jd) {

  diss = 1.49598e8         #Earth-Sun distance (1 AU)
  tmp = moonpos(jd, radian=T)
  ram = tmp$ra
  decm = tmp$dec
  dism = tmp$dis
  tmp = sunpos(jd, radian=T)
  ras = tmp$ra
  decs = tmp$dec
  #cat(ram,decm,dism,ras,decs,sep='\n')
  phi = acos( sin(decs)*sin(decm) + cos(decs)*cos(decm)*cos(ras-ram) )
  inc = atan2( diss * sin(phi), dism - diss*cos(phi) )
  #cat('inc=',inc,'\n')
  k = (1 + cos(inc))/2.
  return(k)
}
