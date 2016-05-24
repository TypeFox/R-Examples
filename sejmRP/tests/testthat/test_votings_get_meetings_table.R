test_that("result of function", {
  expect_equal(class(votings_get_meetings_table("http://www.sejm.gov.pl/Sejm7.nsf/agent.xsp?symbol=posglos&NrKadencji=7")),
    "data.frame")
})

test_that("columns of table", {
  expect_equal(ncol(votings_get_meetings_table("http://www.sejm.gov.pl/Sejm7.nsf/agent.xsp?symbol=posglos&NrKadencji=7")), 3)
})

test_that("rows of table", {
  expect_more_than(nrow(votings_get_meetings_table("http://www.sejm.gov.pl/Sejm8.nsf/agent.xsp?symbol=posglos&NrKadencji=8")), 0)
})