## =============================================================================
##
## Arenstorf orbit
## Standard test problem for nonstiff solvers.
##
## closed trajectory for 3-body problem; two of mass mu and (1-mu)
## and a third body of negligible mass, moving in the same plane
## Hairer et al., 2000
##
## compared with DOPRI.f
##
## =============================================================================

library(deSolve)

#-----------------------------
# the model function
#-----------------------------

Arenstorf <- function(t, y, parms) {
  D1 <- ((y[1] + mu)^2 + y[2]^2)^(3/2)
  D2 <- ((y[1] - (1 - mu))^2 + y[2]^2)^(3/2)
  dy1 <- y[3]
  dy2 <- y[4]
  dy3 <- y[1] + 2*y[4] - (1 - mu)*(y[1] + mu)/D1 - mu*(y[1] - (1 - mu))/D2
  dy4 <- y[2] - 2*y[3] - (1 - mu)*y[2]/D1 - mu*y[2]/D2
  list(c(dy1, dy2, dy3, dy4))
}

#-----------------------------
# parameters, initial values and times
#-----------------------------
mu    <- 0.012277471
yini  <- c(x = 0.994, y = 0, dx = 0, dy = -2.00158510637908252240537862224)
times <- c(seq(from = 0, to = 17, by = 2), 17.0652165601579625588917206249)

#-----------------------------
# solve the model
#-----------------------------
# first for making a graph
system.time({
out <- ode(times = seq(0, 50, 0.1), y = yini, func = Arenstorf, parms = NULL,
  method = rkMethod("ode45"), rtol = 1e-10, atol = 1e-10)
})
plot(out[, c("x", "y")], type = "l", lwd = 2, main = "Arenstorf")

# then for comparison with DOPRI
# (smaller tol than 1e-16 result in numerical problems and very long time)
out <- rk(times = times, y = yini, func = Arenstorf, parms = NULL,
  method = rkMethod("ode45"), rtol = 1e-16, atol = 1e-16)
diagnostics(out)

options(digits = 10)
out[, c("time", "x", "y")]

# this is what DOPRI5 generates with atol=rtol=1e-7:
# X =  0.00    Y =  0.9940000000E+00  0.0000000000E+00    NSTEP =   0
# X =  2.00    Y = -0.5798781411E+00  0.6090775251E+00    NSTEP =  60
# X =  4.00    Y = -0.1983335270E+00  0.1137638086E+01    NSTEP =  73
# X =  6.00    Y = -0.4735743943E+00  0.2239068118E+00    NSTEP =  91
# X =  8.00    Y = -0.1174553350E+01 -0.2759466982E+00    NSTEP = 110
# X = 10.00    Y = -0.8398073466E+00  0.4468302268E+00    NSTEP = 122
# X = 12.00    Y =  0.1314712468E-01 -0.8385751499E+00    NSTEP = 145
# X = 14.00    Y = -0.6031129504E+00 -0.9912598031E+00    NSTEP = 159
# X = 16.00    Y =  0.2427110999E+00 -0.3899948833E+00    NSTEP = 177
# X = XEND     Y =  0.9940021016E+00  0.8911185692E-05
#     tol=0.10D-06   fcn= 1442 step= 240 accpt= 216 rejct= 22


# and this for atol=rtol=1e-17
# X =  0.00    Y =  0.9940000000E+00  0.0000000000E+00    NSTEP =     0
# X =  2.00    Y = -0.5798767232E+00  0.6090783555E+00    NSTEP =  5281
# X =  4.00    Y = -0.1983328832E+00  0.1137637824E+01    NSTEP =  6555
# X =  6.00    Y = -0.4735743108E+00  0.2239077929E+00    NSTEP =  8462
# X =  8.00    Y = -0.1174553507E+01 -0.2759450770E+00    NSTEP = 10272
# X = 10.00    Y = -0.8398071663E+00  0.4468314171E+00    NSTEP = 11505
# X = 12.00    Y =  0.1314377269E-01 -0.8385747019E+00    NSTEP = 13847
# X = 14.00    Y = -0.6031162761E+00 -0.9912585277E+00    NSTEP = 15126
# X = 16.00    Y =  0.2427044376E+00 -0.3899991215E+00    NSTEP = 17184
# X = XEND     Y =  0.9940000000E+00 -0.1966670302E-11
#     tol=0.10D-16   fcn=126836 step=21139 accpt=21137 rejct=    0
