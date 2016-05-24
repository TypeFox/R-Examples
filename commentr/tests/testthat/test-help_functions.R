

context("help_functions")


test_that("full_line", {
  expect_equal(full_line(1), "#") 
  expect_equal(full_line(10), "##########")
  expect_equal(full_line(10, token = "%"), "%%%%%%%%%%")
})


test_that("empty_line", {
  expect_equal(empty_line(1), "# #")
  expect_equal(empty_line(10), "#        #")
  expect_equal(empty_line(10, token = "%"), "%        %")
})


test_that("comment_start", {
  expect_equal(comment_start(), "")
  expect_equal(comment_start(TRUE), "<!--")
})



test_that("comment_end", {
  expect_equal(comment_end(), "")
  expect_equal(comment_end(TRUE), "-->")
})



