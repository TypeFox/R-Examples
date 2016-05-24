page <- read_html("http://www.sejm.gov.pl/Sejm7.nsf/wypowiedz.xsp?posiedzenie=15&dzien=1&wyp=0")
page <- html_nodes(page, ".stenogram")
statements_links <- html_nodes(page, "h2 a")

test_that("result of function", {
  expect_equal(class(statements_get_statements_data(statements_links)), "data.frame")
})

test_that("columns of table", {
  expect_equal(ncol(statements_get_statements_data(statements_links)), 3)
})


test_that("rows of table", {
  expect_more_than(nrow(statements_get_statements_data(statements_links)), 0)
})