grid.int.f <- 0.1
grid.int.psi <- 10

fs <- seq(0, 1, by=grid.int.f)
fs <- c(fs, NA)
psis <- seq(-90, 90, by=grid.int.psi)*pi/180
psis <- c(psis, NA)
gfs   <- outer(fs, psis*0, "+")
gpsis <- outer(fs*0, psis, "+")
gc <- cbind(psi=as.vector(gpsis), f=as.vector(gfs))
gc <- rbind(gc, cbind(psi=as.vector(t(gpsis)), f=as.vector(t(gfs))))

P <- sphere.wedge.to.sphere.cart(gc[,"psi"], gc[,"f"], phi0=135*pi/180)
points3d(P[,"X"], P[,"Y"], P[,"Z"])
lines3d(P[,"X"], P[,"Y"], P[,"Z"])
