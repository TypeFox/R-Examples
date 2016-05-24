## =============================================================================
##
##     This file is adapted from the Test Set for IVP solvers
##     http://www.dm.uniba.it/~testset/
##
##        Andrews'' squeezing mechanism (in index 3 formulation)
##        index 3 DAE of dimension 27
##
## =============================================================================

require(deTestSet)

# -------------------------------------------------------
# problem formulation
# -------------------------------------------------------

# residual function
Andrews <- function (t, y, pars)  {

  with (as.list(pars), {
      sibe <- sin(y[1])
      sith <- sin(y[2])
      siga <- sin(y[3])
      siph <- sin(y[4])
      side <- sin(y[5])
      siom <- sin(y[6])
      siep <- sin(y[7])

      cobe <- cos(y[1])
      coth <- cos(y[2])
      coga <- cos(y[3])
      coph <- cos(y[4])
      code <- cos(y[5])
      coom <- cos(y[6])
      coep <- cos(y[7])

      sibeth <- sin(y[1]+y[2])
      siphde <- sin(y[4]+y[5])
      siomep <- sin(y[6]+y[7])

      cobeth <- cos(y[1]+y[2])
      cophde <- cos(y[4]+y[5])
      coomep <- cos(y[6]+y[7])

      bep <- y[8]
      thp <- y[9]
      php <- y[11]
      dep <- y[12]
      omp <- y[13]
      epp <- y[14]

      m <- matrix(nrow = 7, ncol = 7, data = 0)

      m[1,1] <- m1*ra^2 + m2*(rr^2-2*da*rr*coth+da^2) + i1 + i2
      m[2,2] <- m2*da^2 + i2
      m[3,3] <- m3*(sa^2+sb^2) + i3
      m[4,4] <- m4*(e-ea)^2 + i4
      m[5,5] <- m4*(zt^2+2*zt*(e-ea)*siph+(e-ea)^2) + m5*(ta^2+tb^2) + i4 + i5
      m[6,6] <- m6*(zf-fa)^2 + i6
      m[7,7] <- m6*((zf-fa)^2-2*u*(zf-fa)*siom+u^2) + m7*(ua^2+ub^2) + i6 + i7

      m[2,1] <- m[1,2] <- m2*(da^2-da*rr*coth) + i2
      m[5,4] <- m[4,5] <- m4*((e-ea)^2+zt*(e-ea)*siph) + i4
      m[7,6] <- m[6,7] <- m6*((zf-fa)^2-u*(zf-fa)*siom) + i6

      xd <- sd*coga + sc*siga + xb
      yd <- sd*siga - sc*coga + yb
      lang  <- sqrt ((xd-xc)^2 + (yd-yc)^2)
      force <- - c0 * (lang - l0)/lang
      fx    <- force * (xd-xc)
      fy    <- force * (yd-yc)

      f    <- rep(0,7)
      f[1] <- mom - m2*da*rr*thp*(thp+2*bep)*sith
      f[2] <- m2*da*rr*bep^2*sith
      f[3] <- fx*(sc*coga - sd*siga) + fy*(sd*coga + sc*siga)
      f[4] <- m4*zt*(e-ea)*dep^2*coph
      f[5] <- - m4*zt*(e-ea)*php*(php+2*dep)*coph
      f[6] <- - m6*u*(zf-fa)*epp^2*coom
      f[7] <- m6*u*(zf-fa)*omp*(omp+2*epp)*coom

      gp <- matrix(nrow = 6, ncol = 7, data = 0)

      gp[1,1] <- - rr*sibe + d*sibeth
      gp[1,2] <- d*sibeth
      gp[1,3] <- - ss*coga
      gp[2,1] <- rr*cobe - d*cobeth
      gp[2,2] <- - d*cobeth
      gp[2,3] <- - ss*siga
      gp[3,1] <- - rr*sibe + d*sibeth
      gp[3,2] <- d*sibeth
      gp[3,4] <- - e*cophde
      gp[3,5] <- - e*cophde + zt*side
      gp[4,1] <- rr*cobe - d*cobeth
      gp[4,2] <- - d*cobeth
      gp[4,4] <- - e*siphde
      gp[4,5] <- - e*siphde - zt*code
      gp[5,1] <- - rr*sibe + d*sibeth
      gp[5,2] <- d*sibeth
      gp[5,6] <- zf*siomep
      gp[5,7] <- zf*siomep - u*coep
      gp[6,1] <- rr*cobe - d*cobeth
      gp[6,2] <- - d*cobeth
      gp[6,6] <- - zf*coomep
      gp[6,7] <- - zf*coomep - u*siep

      g    <- rep(0,7)
      g[1] <- rr*cobe - d*cobeth - ss*siga - xb
      g[2] <- rr*sibe - d*sibeth + ss*coga - yb
      g[3] <- rr*cobe - d*cobeth - e*siphde - zt*code - xa
      g[4] <- rr*sibe - d*sibeth + e*cophde - zt*side - ya
      g[5] <- rr*cobe - d*cobeth - zf*coomep - u*siep - xa
      g[6] <- rr*sibe - d*sibeth - zf*siomep + u*coep - ya

      ff        <- rep(0,21)
      ff[1 :14] <- y[8:21]
      ff[15:21] <- - f[1:7] + colSums(m[1:7,1:7 ]* y[15:21]) +
                              colSums(gp[1:6,1:7]* y[22:27])
      ff[22:27] <- g[1:6]

#         ff[1:14]  <- dy[1:14] - ff[1:14]
#         ff[15:27] <- - ff[15:27]
        return(list(ff))
  })
}

# initial conditions
yini <- c(-0.0617138900142764496358948458001, 0,
           0.455279819163070380255912382449,  0.222668390165885884674473185609,
           0.487364979543842550225598953530, -0.222668390165885884674473185609,
           1.23054744454982119249735015568 ,  0,
           0,                                 0,
           0,                                 0,
           0,                                 0,
           14222.4439199541138705911625887,  -10666.8329399655854029433719415,
           0,                                 0,
           0,                                 0,
           0,                                 98.5668703962410896057654982170,
           -6.12268834425566265503114393122,  0,
           0, 0, 0)
           
yprime <- rep(0, times = 27)
yprime[1:14] <- yini[8:21]

# parameters
parameter <- c(m1 = .04325,   m2 = .00365,   m3 = .02373,   m4 = .00706 ,
               m5 = .07050,   m6 = .00706,   m7 = .05498,
               xa = -.06934 , ya = -.00227,
               xb = -0.03635, yb = .03273 ,
               xc = .014 ,    yc = .072,     c0 = 4530,
               i1 = 2.194e-6, i2 = 4.410e-7, i3 = 5.255e-6, i4 = 5.667e-7,
               i5 = 1.169e-5, i6 = 5.667e-7, i7 = 1.912e-5,
               d  = 28e-3,    da = 115e-4,   e = 2e-2,      ea = 1421e-5,
               rr = 7e-3,     ra = 92e-5,    l0 = 7785e-5,  ss = 35e-3,
               sa = 1874e-5,  sb = 1043e-5,  sc = 18e-3,    sd = 2e-2,
               ta = 2308e-5,  tb = 916e-5,
               u = 4e-2,      ua = 1228e-5,  ub = 449e-5,
               zf = 2e-2,     zt = 4e-2,     fa = 1421e-5,   mom = 33e-3)

Mass <- diag (nrow = 27, x = c(rep(1, len = 14), rep (0, len = 13)))


# Check validity of initial conditions
yprime%*%Mass - Andrews(0, yini, parameter)[[1]]


# -------------------------------------------------------
# run at high resolution 
# -------------------------------------------------------
times <- seq(from = 0, to = 0.03, by = 0.001)
ind   <- c(7, 7, 13)
print(system.time(
  MM    <- mebdfi(y = yini, times = times, func = Andrews, mass = Mass,
                  parms = parameter, nind = ind,
                  atol = 1e-10, rtol = 1e-10)
))

print(system.time(
  MM2  <- gamd(y = yini, times = times, func = Andrews, mass = Mass,
                  parms = parameter, nind = ind,
                  atol = 1e-10, rtol = 1e-10)
))

print(system.time(
  MM3    <- radau(y = yini, times = times, func = Andrews, mass = Mass,
                  parms = parameter, nind = ind,
                  atol = 1e-10, rtol = 1e-10)
))

MM[nrow(MM),]
diagnostics(MM)
plot(MM, MM2, MM3)

