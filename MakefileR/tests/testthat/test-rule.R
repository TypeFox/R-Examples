context("rule")

test_that("target-only rule", {
  expect_equal(format(make_rule(".FORCE")), ".FORCE:")
  expect_equal(format(make_rule("target", script = "false")),
               c("target:", "\tfalse"))
})

test_that("target-dep rule", {
  expect_equal(format(make_rule("a", c("b", "c"))), "a: b c")
  expect_equal(format(make_rule("a", "b", c("true", "false"))),
               c("a: b", "\ttrue", "\tfalse"))
  expect_error(make_rule(character()), "target.*required")
})

test_that("appending rules", {
  rules <- list(
    make_rule(".FORCE"),
    make_rule("a", "b")
  )
  expect_equal(makefile(.dots = rules),
               Reduce(c, rules, init = makefile()))
  expect_equal(makefile(.dots = rules),
               makefile() +
                 make_rule(".FORCE") +
                 make_rule("a", "b"))
})

test_that("Printing works as expected", {
  with_mock(
    cat = function(x, sep) x,
    rule <- print(make_rule("a", "b", "true")))
  expect_equal(rule, c("a: b\n", "\ttrue\n"))
})
