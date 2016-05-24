context("comma basics")

test_that('comma works', {
  expect_equal(cc("A", "B"),  "A, B")
  expect_equal(cc("A", c("B", "C")), "A, B, C")
  expect_equal(cc(list("A", "B"), "C"), c("A, B, C"))
})

test_that("data.frame method minimally works", {
  expect_equal(cc(warpbreaks[1, ]), "26, A, L")
  expect_equal(cc(warpbreaks[1:2, ]), "26, A, L, 30, A, L")
  expect_equal(cc(warpbreaks[1:2, "wool"]), "A, A")
})

context("and basics")

test_that("and works", {
  expect_equal(cc_and("A", "B"), "A and B")
  expect_equal(cc_and(c("A", "B"), "C"), "A, B and C")
})

context("or basics")

test_that("or works", {
  expect_equal(cc_or("A", "B"), "A or B")
  expect_equal(cc_or(c("A", "B"), "C"), "A, B or C")
})

context("specials")

test_that('comma special works', {
  expect_equal("A" %,% "B", cc("A","B"))
  expect_equal("A" %,% c("B", "C"), cc("A", "B","C"))
  expect_equal(cc(list("A", "B"), "C"), c("A, B, C"))
})

test_that("and special works", {
  expect_equal("A" %and% "B", cc_and("A","B"))
  expect_equal("A" %and% c("B", "C"), cc_and("A", c("B", "C")))
})

test_that("or special works", {
  expect_equal("A" %or% "B", cc_or("A", "B"))
  expect_equal(c("A", "B") %or% "C", cc_or(c("A", "B"), "C"))
})

context("whitespace")

test_that('whitespace trimming works', {
  expect_equal(cc("  A  ", "  B "), c("A, B"))
})

context("oxford comma")
          
test_that("comma works with cc_or", {
  expect_equal(cc_or("A", "B", "C", oxford = TRUE), "A, B, or C")
  expect_equal(cc_or("A", "B", oxford = TRUE), "A or B")
})

test_that("comma works with cc_and", {
  expect_equal(cc_and("A", "B", oxford = TRUE), "A and B")
  expect_equal(cc_and("A", "B", "C", oxford = TRUE), "A, B, and C")
  expect_equal("A" %and% "B", "A and B")
})
