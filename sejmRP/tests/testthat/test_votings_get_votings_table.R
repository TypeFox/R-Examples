test_that("result of function", {
  expect_equal(class(votings_get_votings_table("http://www.sejm.gov.pl/Sejm7.nsf/agent.xsp?symbol=listaglos&IdDnia=1179")),
    "data.frame")
})

test_that("columns of table", {
  expect_equal(ncol(votings_get_votings_table("http://www.sejm.gov.pl/Sejm7.nsf/agent.xsp?symbol=listaglos&IdDnia=1464")), 3)
})

test_that("rows of table", {
  expect_more_than(nrow(votings_get_votings_table("http://www.sejm.gov.pl/Sejm7.nsf/agent.xsp?symbol=listaglos&IdDnia=1324")), 0)
})