getPDEcsv <-
function(pde = '/home/lees/Site/Santiaguito/pdq.eqs')
  {


eq1 = scan(file=pde, skip=1, sep=",", what=list(yr=0, mo=0, dom=0, t1="", lat=0, lon=0, mag=0, depth=0))

######  parse out the time:
eq1$hr = as.numeric( substr(eq1$t1, 1,2))
eq1$mi = as.numeric( substr(eq1$t1, 3,4))
eq1$sec = as.numeric( substr(eq1$t1, 5,9))

eq1$z = eq1$depth
eq1$jd = getjul(eq1$yr, eq1$mo,  eq1$dom)

invisible(eq1)
  }

