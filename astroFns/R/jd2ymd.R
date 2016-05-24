jd2ymd <-
function(jd) {

# Julian day to year, month, date conversion
# Rounds to nearest Julian date at input
#
# From Fliegel & Van Flandern, Comm. ACM 10, 657 (1968)
#   Note: algorithm uses FORTRAN integer mathematics
#   Also: Explanatory Supplement to the Astronomical Almanac
# A. Harris, U. Maryland Astronomy, 3/17/08

# Convert from 12h to 0h reference, plus a small amount for roundoff
jd <- trunc(jd+0.5)

# Do calculation
L <- jd + 68569
N <- trunc(4 * L/146097)
L <- L - trunc((146097*N + 3)/4)
I <- trunc(4000*(L+1)/1461001)
L <- L - trunc(1461*I/4) + 31
J <- trunc(80*L/2447)
K <- L - trunc(2447*J/80)
L <- trunc(J/11)
J <- J + 2 - 12*L
I <- 100*(N-49) + I + L

ISOdatetime(I, J, K, 0, 0, 0, tz='UTC')

}

