test_that("result of function", {
  expect_equal(class(votes_get_clubs_links("http://www.sejm.gov.pl/Sejm7.nsf/",
    "http://www.sejm.gov.pl/Sejm7.nsf/agent.xsp?symbol=glosowania&NrKadencji=7&NrPosiedzenia=1&NrGlosowania=1")),
    "data.frame")
})

test_that("columns of table", {
  expect_equal(ncol(votes_get_clubs_links("http://www.sejm.gov.pl/Sejm7.nsf/",
    "http://www.sejm.gov.pl/Sejm7.nsf/agent.xsp?symbol=glosowania&NrKadencji=7&NrPosiedzenia=1&NrGlosowania=1")), 2)
})


test_that("rows of table", {
  expect_more_than(nrow(votes_get_clubs_links("http://www.sejm.gov.pl/Sejm7.nsf/",
    "http://www.sejm.gov.pl/Sejm7.nsf/agent.xsp?symbol=glosowania&NrKadencji=7&NrPosiedzenia=1&NrGlosowania=1")), 0)
})