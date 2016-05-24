
context("block_comment")


test_that("Correct output", {
  expect_that(block_comment("hejsan!", width = 20), prints_text("####################"))
  expect_that(block_comment("hejsan!", width = 20), prints_text("#                  #"))
  expect_that(block_comment("hejsan!", width = 20), prints_text("#     hejsan!      #"))
  expect_that(block_comment("hejsan!", html = TRUE), prints_text("<!--"))
  expect_that(block_comment("hejsan!", token = "3"), prints_text("333"))  
  expect_that(block_comment("hejsan!", allign = "left", width = 30), prints_text("# hejsan!                    #"))  
})

