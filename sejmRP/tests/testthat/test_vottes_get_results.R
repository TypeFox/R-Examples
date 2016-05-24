test_that("result of function", {
  expect_equal(class(votes_get_results("http://www.sejm.gov.pl/Sejm7.nsf/agent.xsp?symbol=klubglos&IdGlosowania=37494&KodKlubu=PO")),
    "data.frame")
})

test_that("columns of table", {
  expect_equal(ncol(votes_get_results("http://www.sejm.gov.pl/Sejm7.nsf/agent.xsp?symbol=klubglos&IdGlosowania=42152&KodKlubu=SLD")), 2)
})

test_that("rows of table", {
  expect_more_than(nrow(votes_get_results("http://www.sejm.gov.pl/Sejm7.nsf/agent.xsp?symbol=klubglos&IdGlosowania=32453&KodKlubu=PiS")), 0)
})