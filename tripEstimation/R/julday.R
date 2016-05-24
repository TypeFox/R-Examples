"julday" <-
function(tm) {
  
  yr <- 1900+tm$year
  mon <- tm$mon+1
  day <- tm$mday

  yr[mon<=2] <- yr[mon<=2]-1
  mon[mon<=2] <- mon[mon<=2]+12

  a <- floor(yr/100)
  b <- 2-a+floor(a/4)
  floor(365.25*(yr+4716)+floor(30.6001*(mon+1)))+day+b-1524.5
}

