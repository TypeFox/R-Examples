context("cn")

test_that("minimally works", {
  expect_equal(cn(0, "one", "many"), "one")
  expect_equal(cn(1, "one", "many"), "one")
  expect_equal(cn(1:2, "one", "many"), "many")
})

test_that("%n substitution works", {
  expect_equal(cn("a", "%n number", "%n numbers"), "1 number")
  expect_equal(cn(c("a", "b"), "%n number", "%n numbers"), "2 numbers")
})

test_that("%c substitution works", {
  msg <- c("%n number: %c", "%n numbers: %c")
  expect_equal(cn(7, msg[1], msg[2]), "1 number: 7")
  expect_equal(cn(1:2, msg[1], msg[2]), "2 numbers: 1, 2")
})

test_that("data.frame method works", {
  expect_equal(cn(cars, "one", "many"), "many")
  expect_equal(cn(cars[1, 1, drop = FALSE], "one", "many"), "one")
  expect_equal(cn(cars[, 1, drop = FALSE], "one", "many"), "many")
})

context("cn_and")

test_that("minimally works", {
  expect_equal(cn_and(7, "one", "many"), "one")
  expect_equal(cn_and(1:2, "one", "many"), "many")
})

test_that("%n substitution works", {
  msg <- c("%n number", "%n numbers")
  expect_equal(cn("1", msg[1], msg[2]), "1 number")
  expect_equal(cn(c("1", "2"), msg[1], msg[2]), "2 numbers")
})

test_that("%c substitution works", {
  msg <- c("%n number: %c", "%n numbers: %c")
  expect_equal(cn_and(7, msg[1], msg[2]), "1 number: 7")
  expect_equal(cn_and(1:2, msg[1], msg[2]), "2 numbers: 1 and 2")
})

context("cn_or")

test_that("minimally works", {
  expect_equal(cn_or(7, "one", "many"), "one")
  expect_equal(cn_or(1:2, "one", "many"), "many")
})

test_that("%n substitution works", {
  msg <- c("%n number", "%n numbers")
  expect_equal(cn("1", msg[1], msg[2]), "1 number")
  expect_equal(cn(c("1", "2"), msg[1], msg[2]), "2 numbers")
})

test_that("%c substitution works", {
  msg <- c("%n number: %c", "%n numbers: %c")
  expect_equal(cn_or(7, msg[1], msg[2]), "1 number: 7")
  expect_equal(cn_or(1:2, msg[1], msg[2]), "2 numbers: 1 or 2")
})

test_that("data.frame method minimally works", {
  msg <- c("cars has one row", "cars has many rows")
  expect_equal(cn(cars, msg[1], msg[2]), msg[2])
  expect_equal(cn(cars[1, ], msg[1], msg[2]), msg[1])
  expect_equal(cn(cars[, 1, drop = FALSE], msg[1], msg[2]), msg[2])
})

test_that("%n substitution works with data.frame method", {
  msg <- c("cars has %n row", "cars has %n rows")
  expect_equal(cn(cars[1, ], msg[1], msg[2]), "cars has 1 row")
  expect_equal(cn(cars, msg[1], msg[2]), "cars has 50 rows")
})

test_that("%c substitution works with data.frame method", {
  msg <- c("warpbreaks has %n row: %c", "warpbreaks has %n rows: %c")
  expect_equal(cn(warpbreaks[1, ], msg[1], msg[2]),
               "warpbreaks has 1 row: 26, A, L")
  expect_equal(cn(warpbreaks[1:2, "wool"], msg[1], msg[2]),
               "warpbreaks has 2 rows: A, A")
})

