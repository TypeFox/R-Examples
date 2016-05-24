test_that("result of function", {
  expect_equal(class(votings_get_votings_links("http://www.sejm.gov.pl/Sejm7.nsf/",
    "http://www.sejm.gov.pl/Sejm7.nsf/agent.xsp?symbol=listaglos&IdDnia=1179")),
    "character")
})

test_that("format like link", {
  expect_true(all(stri_detect_regex(votings_get_votings_links("http://www.sejm.gov.pl/Sejm7.nsf/",
    "http://www.sejm.gov.pl/Sejm7.nsf/agent.xsp?symbol=listaglos&IdDnia=1182"),
    "^http")))
})


test_that("length", {
  expect_more_than(length(votings_get_votings_links("http://www.sejm.gov.pl/Sejm8.nsf/",
    "http://www.sejm.gov.pl/Sejm8.nsf/agent.xsp?symbol=listaglos&IdDnia=1491")), 0)
})