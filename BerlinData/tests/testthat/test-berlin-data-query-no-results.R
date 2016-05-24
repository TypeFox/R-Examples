context("berlin_data_query_no_results")

test_that("summary outputs something meaningful", {
  class_obj = structure(list(query = "wat"), 
                        class = "berlin_data_query_no_results")
  expect_that(summary(class_obj), 
              prints_text("Your search did not return any results"))
})