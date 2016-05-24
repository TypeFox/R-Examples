context("Tests list experiment design effect test function")
rm(list=ls())

set.seed(1)

data(affirm)
data(race)

test_that("ict.test works", {
  test.value.affirm <- ict.test(affirm$y, affirm$treat, J = 3, gms = TRUE)
  print(test.value.affirm)
  
  test.value.race <- ict.test(race$y, race$treat, J = 3, gms = TRUE)
  print(test.value.race)
})
