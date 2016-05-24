quadraticPointsOnInterval <- function(t1, t2, n, a) a*seq(0, 1, length=n)^2 + ((t2-t1) - a)*seq(0, 1, length=n) + t1
