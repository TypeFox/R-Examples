
context("Topological sort")

test_that("toplogical sorting works", {

  G <- graph(list(
    "7"  = c("11", "8"),
    "5"  = "11",
    "3"  = c("8", "10"),
    "11" = c("2", "9", "10"),
    "8"  = "9",
    "2"  = character(),
    "9"  = character(),
    "10" = character()
  ))

  topo <- topological_sort(G)

  expect_true(match("7", topo) < match("11", topo))
  expect_true(match("7", topo) < match("8", topo))
  expect_true(match("5", topo) < match("11", topo))
  expect_true(match("3", topo) < match("8", topo))
  expect_true(match("3", topo) < match("10", topo))
  expect_true(match("11", topo) < match("2", topo))
  expect_true(match("11", topo) < match("9", topo))
  expect_true(match("11", topo) < match("10", topo))
  expect_true(match("8", topo) < match("9", topo))

})

test_that("topological sorting works on sparse graphs", {

  G <- graph(list(
    "a" = "b",
    "b" = character(),
    "c" = "d",
    "d" = character(),
    "e" = "f",
    "f" = character()
  ))

  topo <-topological_sort(G)

  expect_true(match("a", topo) < match("b", topo))
  expect_true(match("c", topo) < match("d", topo))
  expect_true(match("e", topo) < match("f", topo))

  topo2 <- topological_sort(
    graph(
      list("a" = "b", "b" = "c", "c" = "d", "d" = "e", "e" = character())
    )
  )

  expect_equal(topo2, c("a", "b", "c", "d", "e"))

  topo3 <- topological_sort(
    graph(list("a" = character(), "b" = character()))
  )

  expect_equal(sort(topo3), c("a", "b"))

})
