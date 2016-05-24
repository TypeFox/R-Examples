# Run but don't test the values. See
# http://cran.r-project.org/doc/manuals/r-release/R-exts.html#Writing-portable-packages
test.fit.density_1 <- function()
{
  model <- WishartModel(100, 400)
  m <- rmatrix(model)

  fitter <- MaximumLikelihoodFit(hint=c(1,1))
  ps <- fit.density(eigen(m), fitter)$par
  #checkEquals(4, ps[1], tolerance=0.1)
  #checkEquals(1, ps[2], tolerance=0.1)
}

test.fit.density_2 <- function()
{
  model <- WishartModel(200, 800, sd=2)
  m <- rmatrix(model)

  fitter <- MaximumLikelihoodFit(hint=c(1,1))
  ps <- fit.density(eigen(cov2cor(m)), fitter)$par
  #checkEquals(4, ps[1], tolerance=0.1)
  #checkEquals(1, ps[2], tolerance=0.1)
}

#test.fit.density_2 <- function()
#{
#  model <- create(WishartModel, 50, 200)
#  m <- rmatrix(model)
#
#  fitter <- create(LegacyFit, hint=c(4,1), kernel='e', adjust=0.2)
#  ps <- fit.density(eigen(cov2cor(m)), fitter)$par
#  checkEquals(4, ps[1], tolerance=0.1)
#  checkEquals(1, ps[2], tolerance=0.1)
#}


# Create a random wishart matrix and verify the cutoff is within the given
# tolerance
# TODO: Test each fit method
test.cutoff_kernel_1 <- function()
{
  model <- WishartModel(50, 200)
  m <- rmatrix(model)
  lp <- cutoff(m)
  #cat("\n")
  #cat("Actual cutoff is",lp,"\n")
  #cat("Theoretical cutoff is",qmp(1, svr=4,var=1),"\n")
  #checkEquals(qmp(1, svr=4, var=1), lp, tolerance=0.1)
}


test.cutoff_kernel_2 <- function()
{
  model <- WishartModel(50, 250)
  m <- rmatrix(model)
  lp <- cutoff(m)
  #cat("\n")
  #cat("Actual cutoff is",lp,"\n")
  #cat("Theoretical cutoff is",qmp(1, svr=5,var=1),"\n")
  #checkEquals(qmp(1, svr=5, var=1), lp, tolerance=0.1)
}


