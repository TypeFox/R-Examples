context("test idig_top_records")

field <- "country"
most <- "united states"
count <- 11
genus <- "acer"
scientificname <- "acer macrophyllum"

test_that("default list of top 10 scientific names returns", {
  testthat::skip_on_cran()
  top <- idig_top_records()
  
  expect_that(top, is_a("list"))
  expect_that(length(top$scientificname), equals(10))
  expect_that(top$itemCount > 20 * 1000 * 1000, is_true())
  
  # Save the number of records in all iDigBio for later tests
  #all_count <- top$itemCount
})

test_that("field and number of tops work", {
  testthat::skip_on_cran()
  top <- idig_top_records(top_fields=c(field), count=count)
  
  expect_that(top, is_a("list"))
  expect_that(length(top[[field]]), equals(count))
  expect_that(top[[field]][[most]][["itemCount"]] > 1000 * 1000, is_true())

  # Deprecating this since Alex changed the erorr behavior to tell you when
  # JSON is bad or field unknown, no longer just spits out all iDigBio
  # Still looking at all of iDigBio, assume things are not changing too fast
  #expect_that(abs(top$itemCount - all_count) < 1000, is_true())
})

test_that("record searches work", {
  testthat::skip_on_cran()
  top <- idig_top_records(rq=list("genus"=genus), top_fields=c(field),
                         count=count)
  
  expect_that(top, is_a("list"))
  expect_that(top$itemCount < 200 * 1000, is_true())
  
  # Save the number of genus records for later tests
  #genus_count <- top$itemCount
})
  
test_that("multiple fields return nested results", {
  testthat::skip_on_cran()
  top <- idig_top_records(rq=list("genus"=genus), top_fields=c(field, 
                         "scientificname"), count=count)
  
  expect_that(top, is_a("list"))
  #expect_that(abs(top$itemCount - genus_count) < 100, is_true())
  expect_that(top[[field]][[most]][["scientificname"]][[scientificname]][["itemCount"]]
              > 1000, is_true())
})