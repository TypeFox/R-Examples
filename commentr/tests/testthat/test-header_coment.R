
context("header_comment")

test_that("Some tests", {
  expect_that(header_comment(), throws_error())
  expect_that(header_comment("ett litet test", author = "Arne Banarne", contact = "no!"), prints_text("ett litet test"))
  expect_that(header_comment("ett litet test", author = "Arne Banarne", contact = "no!"), prints_text("no!"))
  expect_that(header_comment("ett litet test", author = "Arne Banarne", contact = "no!"), prints_text("Arne Banarne"))
  expect_that(header_comment("ett litet test", author = "Arne Banarne", contact = "no!"), prints_text(as.character(Sys.Date())))
})

