grid.coords <- function(grid.int.f=0.5, grid.int.psi=45,
                                              phi0=135*pi/180) {
  ## Funny things happen when f=0
  fs <- seq(1E-6, 1-1E-6, by=grid.int.f - 1E-6)
  psis <- seq(-90, 90, by=grid.int.psi)*pi/180
  
  gfs   <- outer(fs, psis*0, "+")
  gpsis <- outer(fs*0, psis, "+")

  gc <- cbind(psi=as.vector(gpsis), f=as.vector(gfs))
  return(gc)
}

context("Checking wedge coordinates")

test_that("All points lie on sphere", {
  phi0 <- 135*pi/180
  gc <- grid.coords(phi0=phi0)
  P <- sphere.wedge.to.sphere.cart(gc[,"psi"], gc[,"f"], phi0=phi0)
  expect_that(rowSums(P^2), equals(rep(1, nrow(P))))
})

test_that("Some wedge coordinates are converted correctly", {
  phi0 <- 135*pi/180
  expect_that(sphere.wedge.to.sphere.cart(psi=0, f=0.5, phi0=phi0),
              equals(rbind(c(X=0, Y=0, Z=-1))))
  expect_that(sphere.wedge.to.sphere.cart(psi=pi/2, f=0.5, phi0=pi/2),
              equals(rbind(c(X=0, Y=1, Z=0))))
  expect_that(sphere.wedge.to.sphere.cart(psi=-pi/2, f=0.5, phi0=pi/2),
              equals(rbind(c(X=0, Y=-1, Z=0))))
  expect_that(sphere.wedge.to.sphere.cart(psi=-pi/2, f=0, phi0=pi/2),
              equals(rbind(c(X=1, Y=0, Z=0))))
  expect_that(sphere.wedge.to.sphere.cart(psi=0, f=0, phi0=pi/2),
              equals(rbind(c(X=1, Y=0, Z=0))))
  expect_that(sphere.wedge.to.sphere.cart(psi=0, f=1, phi0=pi/2),
              equals(rbind(c(X=-1, Y=0, Z=0))))
})

test_that("Some Cartesian coordinates are converted correctly", {
  expect_that(sphere.cart.to.sphere.wedge(rbind(c(X=0, Y=0, Z=-1)), phi0=pi/2),
              equals(rbind(c(psi=0, f=0.5))))
  expect_that(sphere.cart.to.sphere.wedge(rbind(c(X=0, Y=1, Z=0)), phi0=pi/2),
              equals(rbind(c(psi=pi/2, f=0.5))))
  expect_that(sphere.cart.to.sphere.wedge(cbind(X=-1, Y=0, Z=0), phi0=pi/2)[1,"f"],
              equals(c(f=1)))

})


test_that("Points convert back", {
  phi0 <- 135*pi/180
  gc <- grid.coords(phi0=phi0)
  Pc <- sphere.wedge.to.sphere.cart(gc[,"psi"], gc[,"f"], phi0=phi0)
  Pt <- sphere.cart.to.sphere.wedge(Pc, phi0=phi0)
  expect_that(Pt, equals(gc))
})

context("Checking dualwedge coordinates")
test_that("Some Cartesian coordinates are converted correctly", {
  expect_that(sphere.cart.to.sphere.dualwedge(rbind(c(X=0, Y=0, Z=-1)), phi0=pi/2),
              equals(cbind(fx=0.5, fy=0.5)))
  expect_that(sphere.cart.to.sphere.dualwedge(rbind(c(X=0, Y=1, Z=0)), phi0=pi/2),
              equals(cbind(fx=0.5, fy=0)))
  expect_that(sphere.cart.to.sphere.dualwedge(rbind(c(X=-1, Y=0, Z=0)), phi0=pi/2),
              equals(cbind(fx=1, fy=0.5)))
})
