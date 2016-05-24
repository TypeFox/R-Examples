context("CUSUM charts")

chart <- new("SPCCUSUM",model=SPCModelNormal(Delta=1))


test_that("CUSUM immediate stop",{
    expect_equal(c(spcadjust:::ARL_CUSUM_Markovapprox(1e-8,pobs=pexp)),
                 1,tolerance=1e-3)
})

test_that("getq error messages",{
    expect_error(getq(chart,property="calARL",params=list()),"target")
    expect_error(getq(chart,property="calhitprob",params=list()))
    expect_error(getq(chart,property="calhitprob",params=list(nsteps=100)),"target")

    expect_error(getq(chart,property="calhitprob",params=list(target=0.1)),"nsteps")
    expect_error(getq(chart,property="calhitprob",params=list(target=0.1,nsteps=-1)),"nsteps has to be a positive integer.")
    expect_error(getq(chart,property="calhitprob",params=list(target=0.1,nsteps=0)),"nsteps has to be a positive integer.")
    expect_error(getq(chart,property="calhitprob",params=list(target=0.1,nsteps=1.76)),"nsteps has to be a positive integer.")

    expect_error(getq(chart,property="hitprob",params=list(nsteps=100)),"threshold")
    expect_error(getq(chart,property="ARL",params=list()),"threshold")
    expect_error(getq(chart,property="hitprob",params=list(nsteps=100,threshold=-3)),"Negative threshold.")
    expect_error(getq(chart,property="ARL",params=list(nsteps=100,threshold=-2)),
                 "Negative threshold.")
})


test_that("CUSUM test gridpoints",{
    skip("Extreme case that is not expected to work")
  chart <- new("SPCCUSUM",model=SPCModelNonparCenterScale(Delta=0.1))
  set.seed(347890)
  X <- rexp(100);
  resdefault <- SPCproperty(data=X,nrep=25,chart=chart, property="calARL",params=list(target=10000))
  res150 <- SPCproperty(data=X,nrep=25,chart=chart, property="calARL",params=list(target=10000,gridpoints=150))
  expect_equal(resdefault@res, res150@res, tolerance=0.2)
})

