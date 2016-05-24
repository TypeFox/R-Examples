context("print.HydeNetwork")

data(BlackJack)

test_that("print.HydeNetwork works for full network",
{
  expect_that(print(BlackJack), 
              not(throws_error()))
})

test_that("print.HydeNetwork works for selected nodes",
{
  expect_that(print(BlackJack, dealerFinalPoints, payoff, card5),
              not(throws_error()))
})