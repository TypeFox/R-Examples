test_that("result of function", {
  expect_equal(class(votings_get_date("http://www.sejm.gov.pl/Sejm7.nsf/agent.xsp?symbol=listaglos&IdDnia=1179")),
    "character")
})

test_that("format YYYY-MM-DD", {
  expect_output(votings_get_date("http://www.sejm.gov.pl/Sejm7.nsf/agent.xsp?symbol=listaglos&IdDnia=1179"),
    "[0-9]{4}-[0-9]{2}-[0-9]{2}")
})


test_that("length", {
  expect_equal(length(votings_get_date("http://www.sejm.gov.pl/Sejm8.nsf/agent.xsp?symbol=listaglos&IdDnia=1491")), 1)
})