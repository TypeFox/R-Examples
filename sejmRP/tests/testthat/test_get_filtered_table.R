test_that("result of function", {
  expect_equal(class(get_filtered_votes(host = "services.mini.pw.edu.pl")), "data.frame")
})

test_that("result of function", {
  expect_equal(class(get_filtered_votes(host = "services.mini.pw.edu.pl", clubs = c("PO", "PiS"),
    dates = c("2014-01-01", "2014-12-31"), topics = "referendum")), "data.frame")
})

test_that("result of function", {
  expect_equal(class(get_filtered_votes(host = "services.mini.pw.edu.pl",
    deputies = c("Kopacz Ewa", "Palikot Janusz"))), "data.frame")
})

test_that("columns of table", {
  expect_equal(ncol(get_filtered_votes(host = "services.mini.pw.edu.pl", 
    meetings = c(1, 2))), 9)
})

test_that("columns of table", {
  expect_equal(ncol(get_filtered_votes(host = "services.mini.pw.edu.pl",
    terms_of_office = c(7, 7))), 9)
})

test_that("rows of table", {
  expect_more_than(nrow(get_filtered_votes(host = "services.mini.pw.edu.pl",
    clubs = c("PO", "PiS"), dates = c("2014-01-01", "2014-12-31"),
    topics = "referendum", deputies = c("Kopacz Ewa", "Rostowski"),
    meetings = c(1, 100), votings = c(1, 200))), 0)
})