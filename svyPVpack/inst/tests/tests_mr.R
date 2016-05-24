
# load the survey data. This is the same dataset (as svrepdesign) as in the data folder!
load("svy_example1.RData")

test_that("svyPVpm gives the right output", {

  load("test_svyPVpm.RData")
  expect_that(svyPVpm(by = ~ sex, svydat=svy.exrep, pvs=c("plaus1","plaus2","plaus3")), equals(ex01))
  expect_that(svyPVpm(by = ~ sex + var2, svydat=svy.exrep, pvs=c("plaus1","plaus2","plaus3")), equals(ex02))
  expect_that(svyPVpm(by = ~ sex + var2 + varNA, svydat=svy.exrep, pvs=c("plaus1","plaus2","plaus3")), equals(ex03))
  expect_that(svyPVpm(by = ~ sex + var2 + varNA, svydat=svy.exrep, pvs=c("otherPlaus1","otherPlaus2","otherPlaus3")), equals(ex04))

})


