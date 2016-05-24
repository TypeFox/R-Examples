
context("line_comment")


test_that("Correct output", {
  expect_that(line_comment("Hej pa dig!"), prints_text("Hej pa dig!"))
  expect_that(line_comment("Hejsan svejsan"), prints_text("Hejsan svejsan"))
  expect_that(line_comment("Hejsan svejsan"), prints_text("###########"))  
  expect_that(line_comment("Hejsan svejsan", token = "%"), prints_text("%%%%%%%%"))  
  expect_that(line_comment("Hejsan svejsan", html = TRUE), prints_text("<!--"))  
})
