gal_uvw = function(distance, lsr=F, ra,dec, pmra, pmdec,
  vrad, plx) {
  
  nra = length(ra)
  if(missing(ra) || missing(dec))
    stop('the ra, dec (j2000) position keywords must be supplied (degrees)')

  if(!missing(distance)) {
    if(any(distance<=0)) 
      stop('all distances must be > 0')

    plx = 1e3/distance          #parallax in milli-arcseconds
  }
  else {
    if(missing(plx))
      stop('either a parallax or distance must be specified')
    if(any(plx<=0)) stop('parallaxes must be > 0')
  }
  radeg = 180/pi
  cosd = cos(dec/radeg)
  sind = sin(dec/radeg)
  cosa = cos(ra/radeg)
  sina = sin(ra/radeg)
  k = 4.74047     #equivalent of 1 a.u/yr in km/s   
  a_g = cbind(c( 0.0548755604, +0.4941094279, -0.8676661490),  
         c(0.8734370902, -0.4448296300, -0.1980763734), 
         c(0.4838350155,  0.7469822445, +0.4559837762))
  vec1 = vrad
  vec2 = k*pmra/plx
  vec3 = k*pmdec/plx
  print(( a_g[1,1]*cosa*cosd+a_g[1,2]*sina*cosd+a_g[1,3]*sind)*vec1)
  u = ( a_g[1,1]*cosa*cosd+a_g[1,2]*sina*cosd+a_g[1,3]*sind)*vec1+
      (-a_g[1,1]*sina     +a_g[1,2]*cosa)*vec2+
      (-a_g[1,1]*cosa*sind-a_g[1,2]*sina*sind+a_g[1,3]*cosd)*vec3
  v = ( a_g[2,1]*cosa*cosd+a_g[2,2]*sina*cosd+a_g[2,3]*sind)*vec1+
      (-a_g[2,1]*sina     +a_g[2,2]*cosa)*vec2+
      (-a_g[2,1]*cosa*sind-a_g[2,2]*sina*sind+a_g[2,3]*cosd)*vec3
  w = ( a_g[3,1]*cosa*cosd+a_g[3,2]*sina*cosd+a_g[3,3]*sind)*vec1+
      (-a_g[3,1]*sina     +a_g[3,2]*cosa)*vec2+
      (-a_g[3,1]*cosa*sind-a_g[3,2]*sina*sind+a_g[3,3]*cosd)*vec3

  lsr_vel=c(-8.5,13.38,6.49)
  if(lsr) {
    u = u+lsr_vel[1]
    v = v+lsr_vel[2]
    w = w+lsr_vel[3]
  }

  return(list(u=u,v=v,w=w))
}
