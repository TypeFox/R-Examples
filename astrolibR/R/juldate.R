juldate=function( date) {

  if(length(date)<3)
  stop('illegal date vector - must have a least 3 elements')
  
  date = c(date,rep(0,6-length(date)))

 iy = floor( date[1] ) 
  if(iy<0 )iy = iy +1
  else if(iy==0 ) stop('error - there is no year 0')
  im = floor( date[2] )
 
  day = date[3] + ( date[4] + date[5]/60.0 + date[6]/3600.0) / 24.0
#
if(( im<3 ) ){   #if month is jan or feb, don't include leap day
     iy= iy-1
     im = im+12 
 }
 a = floor(iy/100)
 ry = iy
 jd = floor(ry*0.25) + 365.0*(ry -1860) + floor(30.6001*(im+1.)) + 
      day  - 105.5
if(jd>-100830.5 )jd = jd + 2 - a + floor(a/4)
 
 return(jd)                               
}                                  # juldate
