context("loq_interval performance")

data(ecdata)
data(mfidata)
mfi <- mfidata[mfidata$plate=="plate_1",]
sdf <- data_selection(mfi, ecfile = ecdata)

igmodels <- scluminex("plate_1",sdf$plate_1$standard, sdf$plate_1$background, 
                      lfct="SSl4", bkg="ignore", fmfi="mfi", verbose=FALSE)

consmodels <- scluminex("plate_1",sdf$plate_1$standard, sdf$plate_1$background, 
                      lfct="SSl4", bkg="constraint", fmfi="mfi", verbose=FALSE)



test_that("no input error", {
  expect_that(loq_interval(NA), throws_error())
})

test_that("bad arguments combination ignore model", {
  expect_that(loq_interval(igmodels, low.asymp=2, lowci=2), throws_error())
  expect_that(loq_interval(igmodels, low.asymp=2, lowci=Inf), throws_error())
  expect_that(loq_interval(igmodels, high.asymp=2, highci=-Inf), throws_error())
  expect_that(loq_interval(igmodels, high.asymp=2, highci=-Inf), throws_error())
  expect_that(loq_interval(igmodels, low.asymp=NULL, lowci=Inf), throws_error())
  expect_that(loq_interval(igmodels, high.asymp=NULL, highci=-Inf), throws_error())
  expect_that(loq_interval(igmodels, high.asymp=40), throws_error())
  expect_that(loq_interval(igmodels, low.asymp=40), throws_error())
})

test_that("bad arguments combination constraint model", {
  expect_that(loq_interval(consmodels, low.asymp=2, lowci=2), throws_error())
  expect_that(loq_interval(consmodels, low.asymp=2, lowci=Inf), throws_error())
  expect_that(loq_interval(consmodels, high.asymp=2, highci=-Inf), throws_error())
  expect_that(loq_interval(consmodels, high.asymp=2, highci=-Inf), throws_error())
  expect_that(loq_interval(consmodels, low.asymp=NULL, lowci=Inf), throws_error())
  expect_that(loq_interval(consmodels, high.asymp=NULL, highci=-Inf), throws_error())
})

test_that("check arguments specific constraint model", {
  expect_that(loq_interval(consmodels, low.asymp=2), is_a("loq"))
  expect_that(loq_interval(consmodels, low.asymp=2, high.asymp=8), throws_error())
  expect_that(loq_interval(consmodels, high.asymp=40), throws_error())
  expect_that(loq_interval(consmodels, low.asymp=40), is_a("loq"))
})





