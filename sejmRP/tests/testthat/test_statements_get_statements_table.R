test_that("result of function", {
  expect_equal(class(statements_get_statements_table("http://www.sejm.gov.pl/Sejm7.nsf/posiedzenie.xsp?posiedzenie=99&dzien=2")),
    "data.frame")
})

test_that("columns of table", {
  expect_equal(ncol(statements_get_statements_table("http://www.sejm.gov.pl/Sejm7.nsf/posiedzenie.xsp?posiedzenie=50&dzien=1")), 3)
})


test_that("rows of table", {
  expect_more_than(nrow(statements_get_statements_table("http://www.sejm.gov.pl/Sejm8.nsf/posiedzenie.xsp?view=1&posiedzenie=2&dzien=1")), 0)
})