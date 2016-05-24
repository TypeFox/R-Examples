co_refract=function(
  a, altitude=0, pressure,
  temperature, to_observed, epsilon=0.25) {
                                        #epsilon=0.25: accuracy of iteration for observed=1 case, in arcseconds
  alpha = 0.0065 # temp lapse rate [deg c per meter]

  if(missing(pressure)) 
    pressure = 1010.*(1-6.5/288000*altitude)^5.255
	
  if (missing(temperature)) {
    if(altitude>11000)
      temperature = 211.5
    else
      temperature = 283.0 - alpha*altitude
  }

  if(missing(to_observed)){
    aout = a - co_refract_forward(a,p=pressure,t=temperature)
  }
  else {
    aout = a*0.
    na = length(a)
    p = pressure + a*0.
    t = temperature + a*0.
    for(i in 1:na) {
                                        #calculate initial refraction guess
      dr = co_refract_forward(a[i],p=p[i],t=t[i])
      cur = a[i] + dr # guess of observed location
      while(t) {
        last = cur
        dr = co_refract_forward(cur,p=p[i],t=t[i])
        cur= a[i] + dr
        if(abs(last-cur)*3600.<epsilon) break
      }
      aout[i] = cur
    }
  }
  return (aout)
}
