test_that("result of function", {
  expect_equal(class(statements_get_statement(
    "http://www.sejm.gov.pl/Sejm7.nsf/wypowiedz.xsp?posiedzenie=15&dzien=1&wyp=008")),
    "character")
})

test_that("length of vector", {
  expect_equal(length(statements_get_statement(
    "http://www.sejm.gov.pl/Sejm8.nsf/wypowiedz.xsp?posiedzenie=2&dzien=1&wyp=021")), 1)
})