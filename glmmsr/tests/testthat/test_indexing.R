library(glmmsr)
context("Indexing")

test_that("checks dimension of indexing correctly", {
  subexpr1 <- quote(ability[player1] - ability[player2])
  subexpr2 <- quote(ability[player1, player2] - ability[player2, player1])
  subexpr3 <- quote(ability[player1] + stuff[m1, m2, m3])
  subexpr4 <- quote(ability[player1] + ability[player1, player2])
  expect_equal(find_dim_sub(subexpr1, "ability"), 1L)
  expect_equal(find_dim_sub(subexpr2, "ability"), 2L)
  expect_equal(find_dim_sub(subexpr3, "ability"), 1L)
  expect_equal(find_dim_sub(subexpr3, "stuff"), 3L)
  expect_error(find_dim_sub(subexpr4, "ability"),
               "indexed inconsistently")
})


test_that("deduces indices for subform correctly", {
  sub1 <- list(subvar = "ability",
               subform = formula(ability[player] ~ 0 + (1 | player)),
               subexpr = quote(ability[player1] - ability[player2]))
  player1 <- factor(c("a", "b", "c"), levels = c("a", "b", "c", "d"))
  player2 <- factor(c("b", "c", "d"), levels = c("a", "b", "c", "d"))
  data1 <- list(player1 = player1, player2 = player2)
  indices1 <- find_indices_subform(sub1, data1)$indices_subform
  expect_equal(names(indices1), "player")
  expect_equal(length(indices1$player), 4L)

  sub2 <- list(subvar = "ability",
               subform = formula(ability[p, m] ~ 0 + x[m] + (1 | p)),
               subexpr = quote(ability[player1, match]
                               - ability[player2, match]))
  data2 <- list(player1 = player1,
                player2 = player2,
                match = c(1, 2, 3))
  indices2 <- find_indices_subform(sub2, data2)$indices_subform
  expect_equal(names(indices2), c("p", "m"))
  expect_equal(length(indices2$p), 4L)
  expect_equal(length(indices2$m), 3L)

  player1 <- factor(c("a", "b", "c"), levels = c("a", "b", "c"))
  player2 <- factor(c("b", "c", "d"), levels = c("b", "c", "d"))
  data3 <- list(player1 = player1, player2 = player2)
  expect_error(find_indices_subform(sub1, data3)$indices_subform,
               "Indexing factors must have identical levels")

})

test_that("converts from multi-indexing to single indexing correctly", {
  indices <- list("a" = c(1, 2, 3), "b" = c(1, 2), "c" = c(1, 2, 3, 4))
  indices_exp <- expand.grid(indices)
  ids <- multi_to_flat(t(indices_exp[c(3, 5, 11, 23),]), indices)
  expect_equal(ids, c(3, 5, 11, 23))
})
