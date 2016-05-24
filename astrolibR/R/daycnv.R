daycnv = function( xjd) {

 jd = floor(xjd)                         #Truncate to integral day
 frac = xjd - jd + 0.5          #Fractional part of calendar day
 after_noon =which(frac>=1.0)

 frac[after_noon] = frac[after_noon] - 1.0
 jd[after_noon] = jd[after_noon] + 1

 hr = frac*24.0
 L = jd + 68569
 n = 4*L / 146097
 L = L - (146097*n + 3) / 4
 yr = 4000*(L+1) / 1461001
 L = L - 1461*yr / 4 + 31        #1461 = 365.25 * 4
 mn = 80*L / 2447
 day = L - 2447*mn / 80
 L = mn/11
 mn = mn + 2 - 12*L
 yr = 100*(n-49) + yr + L
 return (list(yr=round(yr),mn=round(mn),day=round(day),hr=hr))
}
