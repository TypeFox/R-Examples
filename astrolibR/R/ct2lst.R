#source('jdcnv.R')

ct2lst=function(lng, tz, tme, day, mon, year) {

 if(missing(day) || missing(mon) || missing(year)) {
   jd = tme;
 }
 else {
   time = tme - tz
   jd = jdcnv(year, mon, day, time) 
 }

 c1 = c(280.46061837e0,
       360.98564736629e0,
       0.000387933e0,
       38710000.0)
 jd2000 = 2451545.0
 t0 = jd - jd2000
 t = t0/36525

 theta = c1[1] + (c1[2] * t0) + t^2*(c1[3] - t/ c1[4] )

 lst = ( theta + lng)/15
 neg = (lst<0);
 lst[neg] = 24 + (lst[neg] %% 24)
 lst = lst %% 24
 return(lst)
}

