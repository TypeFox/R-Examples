library(glmmsr)
context("Modifying subexpressions")

test_that("modifies subexpr for model frame indexing correctly", {
  subexpr1 <- quote(ability[player1] - ability[player2])
  subexpr2 <- quote(ability[player1] + stuff[id])
  expect_equal(modify_subexpr(subexpr1, "ability"),
               quote(`-fr`(`[fr`(ability, player1), `[fr`(ability, player2))))
  expect_equal(modify_subexpr(subexpr2, "ability"),
               quote(`+fr`(`[fr`(ability, player1), stuff[id])))
  subexpr3 <- quote(ability[player1, match] - ability[player2, match])
  subexpr3_mod <- modify_subexpr(subexpr3, "ability")
  to_flatten <- extract_to_flatten(subexpr3, "ability")
  to_flatten_manual <- list(c("player1", "match"), c("player2", "match"))
  expect_equal(to_flatten, to_flatten_manual)
  subexpr4 <- quote(ability[player1, match, stuff])
  to_flatten_4 <- extract_to_flatten(subexpr4, "ability")
  expect_equal(to_flatten_4, list(c("player1", "match", "stuff")))
  subexpr5 <- quote(ability[player1, match, stuff] + other[one, two, three])
  to_flatten_5 <- extract_to_flatten(subexpr5, "ability")
  expect_equal(to_flatten_5, list(c("player1", "match", "stuff")))
})
