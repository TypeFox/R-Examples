
###########
context("Ace Estimation")
###########
test_that("CreateAceEstimate -Plain", {
  aSquared <- .5
  cSquared  <- .3
  eSquared <- .2
  caseCount <- 20
  
  est <- CreateAceEstimate(aSquared, cSquared, eSquared, caseCount)
  expect_equal(object=slot(est, "ASquared"), expected=aSquared, scale=1)
  expect_equal(object=slot(est, "CSquared"), expected=cSquared, scale=1)
  expect_equal(object=slot(est, "ESquared"), expected=eSquared, scale=1)
  expect_equal(object=slot(est, "CaseCount"), expected=caseCount, scale=1)
  expect_true(object=slot(est, "Unity"))
  expect_true(object=slot(est, "WithinBounds"))
})
test_that("CreateAceEstimate -Not Unity", {
  aSquared <- .5
  cSquared  <- .2
  eSquared <- .1
  caseCount <- 20
  
  est <- CreateAceEstimate(aSquared, cSquared, eSquared, caseCount)
  expect_equal(object=slot(est, "ASquared"), expected=aSquared, scale=1)
  expect_equal(object=slot(est, "CSquared"), expected=cSquared, scale=1)
  expect_equal(object=slot(est, "ESquared"), expected=eSquared, scale=1)
  expect_equal(object=slot(est, "CaseCount"), expected=caseCount, scale=1)
  expect_false(object=slot(est, "Unity"))
  expect_true(object=slot(est, "WithinBounds"))
})
test_that("CreateAceEstimate -Outside bounds", {
  aSquared <- .7
  cSquared  <- .5
  eSquared <- -.2
  caseCount <- 20
  
  est <- CreateAceEstimate(aSquared, cSquared, eSquared, caseCount)
  expect_equal(object=slot(est, "ASquared"), expected=aSquared, scale=1)
  expect_equal(object=slot(est, "CSquared"), expected=cSquared, scale=1)
  expect_equal(object=slot(est, "ESquared"), expected=eSquared, scale=1)
  expect_equal(object=slot(est, "CaseCount"), expected=caseCount, scale=1)
  expect_true(object=slot(est, "Unity"))
  expect_false(object=slot(est, "WithinBounds"))
})
test_that("CreateAceEstimate -Outside bounds & no unity", {
  aSquared <- .4
  cSquared  <- .5
  eSquared <- -.2
  caseCount <- 20
  
  est <- CreateAceEstimate(aSquared, cSquared, eSquared, caseCount)
  expect_equal(object=slot(est, "ASquared"), expected=aSquared, scale=1)
  expect_equal(object=slot(est, "CSquared"), expected=cSquared, scale=1)
  expect_equal(object=slot(est, "ESquared"), expected=eSquared, scale=1)
  expect_equal(object=slot(est, "CaseCount"), expected=caseCount, scale=1)
  expect_false(object=slot(est, "Unity"))
  expect_false(object=slot(est, "WithinBounds"))
})

test_that("print AceEstimate shows something", {
  est <- CreateAceEstimate(aSquared=.5, cSquared= .2, eSquared= .4, caseCount=4)
  expect_that(print(est), prints_text(regexp="[:alnum:]"))
})

# 
# aSquared <- .5
# cSquared  <- .2
# eSquared <- .1
# caseCount <- 20
# componentSum <- aSquared + cSquared + eSquared
# unity <- ( abs(componentSum - 1.0) < 0 )
# withinBounds <- (0 <= min(aSquared, cSquared, eSquared)) && (max( aSquared, cSquared, eSquared) <= 1)
# # est <-new("AceEstimate", aSquared, cSquared, eSquared, caseCount, unity, withinBounds) 
# # est@ASquared
# # est
# # show(est)
# # print(est)
# # 
# est2 <- CreateAceEstimate(.5, .2, .1, 4)
# # est2@ASquared
# # est2
# # show(est2)
# # # showClass("AceEstimate")
# 
# expect_that(print(est2), prints_text(regexp="[:alnum:]"))
#           
# expect_that(message(print(est2)), shows_message())
# expect_output(print(est2))
# expect_that(print(est2), prints_text())
# expect_that(print(est2), prints_text("\"Results of ACE estimation: [show]\"\nASquared  CSquared  ESquared CaseCount\n0.5       0.2       0.1       4.0)"))

# showMethods(GetEstimate)
# showMethods(print)
# showMethods(show)
