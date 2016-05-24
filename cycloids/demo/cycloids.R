## Demo for Package cycloids
## Author: Peter Biber, November 2013

## A classic spirograph pattern
## Varying parameters (here: lambda) within a loop often gives
## nice results.

op <- par(mar = c(0,0,0,0), bg = "black") # no plot margins
plot.new()
plot.window(asp = 1, xlim = c(-4.5, 4.5), ylim = c(-4.5, 4.5))
lambdax <- seq(0.85, by = -0.05, length.out = 14)
ccol <- rep(c("steelblue", "steelblue", "red", "red"), 4)
# draw fourteen hypotrochoids with decreasing lambda
for (i in c(1:14)) {
     z <- zykloid(5, 3, lambdax[i])
     lines(y ~ x, data = z, type = "l", col = ccol[i])
} # for i
par(op) # set graphics parameters back to original values


## This pattern keeps the outer radius of the cycloid constant
## while modifying lambda.

op <- par(mar = c(0,0,0,0), bg = "black")
plot.new()
plot.window(asp = 1, xlim = c(-1, 1), ylim = c(-1, 1))
lambdax <- seq(2, 0.0, -0.05) # Note: some lambdas are greater than 1
ccol <- rep(c("lightblue", "lightblue", "yellow", "yellow", "yellow"), 9)
for (ll in c(1:length(lambdax))) {
     z <- zykloid.scaleP(A = 7, a = 5, hypo = TRUE, lambda = lambdax[ll])
     lines(y ~ x, data = z, col = ccol[ll])
} # for ll
par(op) # set graphics parameters back to original values


## Spiky Flower with zykloid.scaleA and zykloid

op <- par(mar = c(0,0,0,0), bg = "black")
plot.new()
plot.window(asp = 1, xlim = c(-150, 150), ylim = c(-150, 150))
z <- zykloid.scaleA(A = 90, a = 32, lambda = 1, Radius = 150, hypo = TRUE)
lines(y ~ x, data = z, col = "lightblue")
for (ll in seq(2, 0.8, -0.4)) {
     if (ll == 2) ccol <- "royalblue"
     else         ccol <- "plum"
     z <- zykloid(A = 90, a = 32, lambda = ll, hypo = TRUE, steps = 360, start = pi/2)
     lines(y ~ x, data = z, col = ccol)
} # for ll
par(op)


## Looks like a passion flower

op <- par(mar = c(0,0,0,0), bg = "black")
plot.new()
plot.window(asp = 1, xlim = c(-23, 23), ylim = c(-23, 23))
ll   <- seq(2, 0, -0.2)
ccol <- rep(c("lightblue", "lightgreen", "yellow", "yellow",
              "yellow"), 2)
for (i in c(1:length(ll))) {
     z <- zykloid(A = 15, a = 7, lambda = ll[i], hypo = TRUE)
     lines(y ~ x, data = z, col = ccol[i])
} # for i
par(op)


## Dense hypotrochoids

op <- par(mar = c(0,0,0,0), bg = "black")
plot.new()
plot.window(asp = 1, xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5))
m <- zykloid(A = 90, a = 89, lambda = 0.01)
lines(y ~ x, data = m, col = "grey")
m <- zykloid(A = 90, a = 89, lambda = 0.02)
lines(y ~ x, data = m, col = "red")
m <- zykloid(A = 90, a = 89, lambda = 0.015)
lines(y ~ x, data = m, col = "blue")
par(op)


## Fragile star

op <- par(mar = c(0,0,0,0), bg = "black")
plot.new()
plot.window(asp = 1, xlim = c(-14, 14), ylim = c(-14, 14))
l.max <- 1.6
l.min <- 0.1
ll <- seq(l.max, l.min, by = -1 * (l.max - l.min)/30)
n  <- length(ll)
ccol <- rainbow(n, start = 2/3, end = 1)
for (i in c(1:n)) {
    m <- zykloid(A = 9, a = 8, lambda = ll[i])
    lines(y ~ x, data = m, type = "l", col = ccol[i])
}  # for i
par(op)


## In this example, RadiusA depends on the cosine of the x-coordinate
## of the fixed circle's centre

op <- par(mar = c(0,0,0,0), bg = "black")
plot.new()
plot.window(asp = 1, xlim = c(-8, 8), ylim = c(-0.5, 0.5))
ctrx <- seq(-2*pi, 2*pi, pi/10)
ccol <- rainbow(length(ctrx))
for(i in c(1:length(ctrx))) {
    zzz <- zykloid.scaleA(A = 9, a = 7, hypo = TRUE, Cx = ctrx[i],
                          Cy = -ctrx[i], lambda = 0.9,
                          RadiusA = 1.5 + cos(ctrx[i]), start = -pi/4)
    lines(y ~ x, data = zzz, col = ccol[i])
} # for i
par(op)



## Geometric degression of RadiusA makes a nice star

op <- par(mar = c(0,0,0,0), bg = "black")
plot.new()
plot.window(asp = 1, xlim = c(-10, 10), ylim = c(-10, 10))
rad <- 10
n <- 60
ccol <- heat.colors(n)
for(i in c(1:n)) {
    if (i/2 != floor(i/2)) { sstart = pi/2 }
    else                   { sstart = pi/4 }
    zzz <- zykloid.scaleA(A = 4, a = 3, RadiusA = rad, lambda = 1,
                          start = sstart)
    lines(y ~ x, data = zzz, col = ccol[i])
    rad <- rad * 0.9
} # for i
par(op)



## A windmill

op <- par(mar = c(0,0,0,0), bg = "black")
plot.new()
plot.window(asp = 1, xlim = c(-1.4, 1.4), ylim = c(-1.4, 1.4))
rrad <- sqrt(seq(0.1, 2, 0.1))
n    <- length(rrad)
ccol <- rainbow(n, start = 0, end = 0.3)
for(i in c(1:n)) {
    zzz <- zykloid.scaleA(A = 7, a = 3, RadiusA = rrad[i],
           hypo = TRUE, lambda = 1.1,
           start = pi/2 - (1*pi/7 - (i - 1) * 2*pi/(7 * n)))
    lines(y ~ x, data = zzz, col = ccol[n + 1 - i])
} # for i
par(op)



## Advanced Example: A series of cycloids with their centres
## located on a logarithmic spiral

op <- par(mar = c(0,0,0,0), bg = "black")
plot.new()
plot.window(asp = 1, xlim = c(-50, 50), ylim = c(-50, 50))
a     <- 1/32     # spiral's scaling constant
alpha <- pi/20    # spiral's slope angle
sphi  <- seq(0, 18 * pi, pi/25)   # series of angles for cycloid centres
rad  <- a * exp(tan(alpha)*sphi)  # corresponding spiral radii
spx  <- rad * cos(sphi)           # corresponding x-coordinates
spy  <- rad *sin(sphi)            # corresponding y-coordinates
n    <- length(sphi)
ccol <- rainbow(n, start = 2/3, end = 1/2)
for (i in c(1:n)) {
     czc <- zykloid.scaleA(A = 3, a = 1, lambda = 1.5,
            Cx = spx[i], Cy = spy[i],
            RadiusA = rad[i]/2.5, # cycloid radii depends on spiral radii
            start = pi + sphi[i]) # angle cycloid towards spiral centre
     lines(y ~ x, data = czc, col = ccol[i])
} # for i
par(op)


## Pentagram by constructing a hypocycloid and an epicycloid
## with the same outer radius and scaling this radius exponentially

op <- par(mar = c(0,0,0,0), bg = "black")
plot.new()
plot.window(asp = 1, xlim = c(-40, 40), ylim = c(-40, 40))
n <- 20
ccol <- heat.colors(n)
for(i in c(1:n)) {
    zzz <- zykloid.scaleAa(A = 5, a = 2,
           RadiusAa = 38*exp(-0.05*(i-1)), hypo = FALSE, lambda = 1)
    lines(y ~ x, data = zzz, col = ccol[i])
    zzz <- zykloid.scaleAa(A = 5, a = 2,
           RadiusAa = 38*exp(-0.05*(i-1)), hypo = TRUE,  lambda = 1)
    lines(y ~ x, data = zzz, col = ccol[i])
} # for i
par(op)



## Psychedelic star by modifying lambda while keeping the outer
## radius constant

op <- par(mar = c(0,0,0,0), bg = "black")
plot.new()
plot.window(asp = 1, xlim = c(-5, 5), ylim = c(-5, 5))
llam <- seq(0, 8, 0.2)
ccol <- terrain.colors(length(llam))
for(i in c(1:length(llam))) {
    zzz <- zykloid.scaleAa(A = 5, a = 1, RadiusAa = 4.5,
           hypo = FALSE, lambda = llam[i])
    lines(y ~ x, data = zzz, col = ccol[i])
} # for i
par(op)


## Cool Disk by scaling the start angle with an
## exponential function ...

op <- par(mar = c(0,0,0,0), bg = "black")
plot.new()
plot.window(asp = 1, xlim = c(-11, 11), ylim = c(-11, 11))
n <- 30
ccol <- topo.colors(n)
for(i in c(1:n)) {
    zzz <- zykloid.scaleP(A = 3, a = 1, RadiusP = 6, lambda = 1,
           start = 2*pi/3 * exp(-0.1 * (i - 1)), hypo = FALSE)
    lines(y ~ x, data = zzz, col = ccol[i])
} # for i
par(op)



## ... the free space in the centre could be filled with
## the corresponding hypocycloid ...

op <- par(mar = c(0,0,0,0), bg = "black")
plot.new()
plot.window(asp = 1, xlim = c(-11, 11), ylim = c(-11, 11))
n <- 30
ccol <- topo.colors(n)
for(i in c(1:n)) {
    zzz <- zykloid.scaleP(A = 3, a = 1, RadiusP = 6, lambda = 1,
           start = 2*pi/3 * exp(-0.1 * (i - 1)), hypo = FALSE)
    lines(y ~ x, data = zzz, col = ccol[i])
    zzz <- zykloid.scaleP(A = 3, a = 1, RadiusP = 6, lambda = 1,
           start = 2*pi/3 * exp(-0.1 * (i - 1)), hypo = TRUE)
    lines(y ~ x, data = zzz, col = ccol[i])
} # for i
par(op)



## ... or the same ring again and again.

op <- par(mar = c(0,0,0,0), bg = "black")
plot.new()
plot.window(asp = 1, xlim = c(-11, 11), ylim = c(-11, 11))
n <- 30
ccol <- topo.colors(n)
rad <- 6
for(g in c(1:7)) {
    for(i in c(1:n)) {
        zzz <- zykloid.scaleP(A = 3, a = 1, RadiusP = rad,
               lambda = 1, start = 2*pi/3 * exp(-0.1 * (i - 1)),
               hypo = FALSE)
        lines(y ~ x, data = zzz, col = ccol[i])
    } # for i
    rad <- rad * 3/5
} # for g
par(op)



## Cauliflower pattern. Here, an exponential function is used
## for scaling the radius of the circle the cycloid's loops
## are on.

op <- par(mar = c(0,0,0,0), bg = "black")
plot.new()
plot.window(asp = 1, xlim = c(-22, 22), ylim = c(-22, 22))
n <- 15
dcol <- heat.colors(n)
for(i in c(1:n)) {
    lambdax <- seq(2.0, 2.2, 0.1)
    for(j in c(1:length(lambdax))) {
        zzz <- zykloid.scaleP(A = 11, a = 1,
               RadiusP = 15 * exp(-0.3 * (i - 1)),
               lambda = lambdax[j], hypo = FALSE,
               start = pi/2 + (i - 1)*pi/11)
        if(j/2 == floor(j/2)) { colx <- "blue" }
        else                  { colx <- dcol[n + 1 - i] }
        lines(y ~ x, data = zzz, col = colx)
    } # for j
} # for i
par(op)



## Sparkling star

op <- par(mar = c(0,0,0,0), bg = "black")
plot.new()
plot.window(asp = 1, xlim = c(-15, 15), ylim = c(-15, 15))
llam <- seq(0, 8, 0.2)
ccol <- rainbow(length(llam), start = 2/3, end = 1/3)
for(i in c(1:length(llam))) {
    zzz <- zykloid.scaleP(A = 5, a = 1, RadiusP = 2.1,
           hypo = FALSE, lambda = llam[i], start = pi/5)
    lines(y ~ x, data = zzz, col = ccol[i])
} # for i
par(op)






