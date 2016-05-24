test_that("result of function", {
  expect_equal(class(votes_match_deputies_ids("sejmrp", "reader", "qux94874", "services.mini.pw.edu.pl",
    "http://www.sejm.gov.pl/Sejm7.nsf/agent.xsp?symbol=klubglos&IdGlosowania=37494&KodKlubu=PO", 7,
    .Platform$OS.type == "windows")),
    "data.frame")
})

test_that("columns of table", {
  expect_equal(ncol(votes_match_deputies_ids("sejmrp", "reader", "qux94874", "services.mini.pw.edu.pl",
    "http://www.sejm.gov.pl/Sejm7.nsf/agent.xsp?symbol=klubglos&IdGlosowania=37494&KodKlubu=SLD", 7,
    .Platform$OS.type == "windows")), 3)
})

test_that("rows of table", {
  expect_more_than(nrow(votes_match_deputies_ids("sejmrp", "reader", "qux94874", "services.mini.pw.edu.pl",
    "http://www.sejm.gov.pl/Sejm7.nsf/agent.xsp?symbol=klubglos&IdGlosowania=37494&KodKlubu=PiS", 7,
    .Platform$OS.type == "windows")), 0)
})