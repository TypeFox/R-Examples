#source('ten.R')
#source('nutate.R')
#source('sunpos.R')

co_aberration = function( jd, ra, dec, eps) {
  
  d2r = pi/180
  T = (jd -2451545.0)/36525.0 # julian centuries from J2000 of jd.
  if(missing(eps)){ # must calculate obliquity of ecliptic
    njd = length(jd)
    d_psi = numeric(njd)
    d_epsilon = d_psi
    for(i in 1:njd) {
      tmp=nutate(jd[i])
      dp=tmp$nut_long
      de=tmp$nut_obliq # d_psi and d_epsilon in degrees
      d_psi[i] = dp
      d_epsilon[i] = de
    }
    eps0 = ten(23,26,21.448)*3600 - 46.8150*T - 0.00059*T^2 +  
      0.001813*T^3
    eps = (eps0 + d_epsilon)/3600.*d2r # true obliquity of the ecliptic 
  }
  
  tmp = sunpos(jd)
  sunra = tmp$ra
  sundec = tmp$dec
  sunlon = tmp$longmed
  
  e = 0.016708634 - 0.000042037*T - 0.0000001267*T^2
  pi = 102.93735 + 1.71946*T + 0.00046*T^2 
  k = 20.49552 #constant of aberration, in arcseconds
  cd = cos(dec*d2r) ; sd = sin(dec*d2r)
  ce = cos(eps) ; te = tan(eps)
  cp = cos(pi*d2r) ; sp = sin(pi*d2r)
  cs = cos(sunlon*d2r) ; ss = sin(sunlon*d2r)
  ca = cos(ra*d2r) ; sa = sin(ra*d2r)
  term1 = (ca*cs*ce+sa*ss)/cd
  term2 = (ca*cp*ce+sa*sp)/cd
  term3 = (cs*ce*(te*cd-sa*sd)+ca*sd*ss)
  term4 = (cp*ce*(te*cd-sa*sd)+ca*sd*sp)
  d_ra = -k * term1 + e*k * term2
  d_dec = -k * term3 + e*k * term4

  return(list(d_ra=d_ra,d_dec=d_dec,eps=eps))
}

