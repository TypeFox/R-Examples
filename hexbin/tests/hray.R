library(hexbin)

set.seed(572)

x <- rnorm(100)
y <- rnorm(100)
val <- rnorm(100)
inc <- abs(rnorm(100,sd = .3))
loB <- val-inc
hiB <- val+inc

if(exists("hray", mode="function")) { # 'real soon now'

## no confidence bounds
plot(x,y,type = 'n')
hray(x,y,val)

## confidence bounds
plot(x,y,type = 'n')
hray(x,y,val, lo = loB, hi = hiB)

## clockwise orientation
plot(x,y,type = 'n')
hray(x,y,val, loB, hiB, clockwise = TRUE)

## no tics and small filled dots
plot(x,y,type = 'n')
hray(x,y,val, loB, hiB, ticlength = FALSE,
     dotside = 20, dotlength = .025, dotden = -1)

}
