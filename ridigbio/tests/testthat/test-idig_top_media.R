context("test idig_top_media")

field <- "recordset"
most <- "7450a9e3-ef95-4f9e-8260-09b498d2c5e6"
count <- 11
dqs <- 0
scientificname <- "acer macrophyllum"

test_that("default list of top 10 scientific names returns", {
  testthat::skip_on_cran()
  top <- idig_top_media()
  
  expect_that(top, is_a("list"))
  #FIXME: Alex should be changing the default summary to hasSpecimen soon
  #expect_that(length(top[["hasSpecimen"]]), equals(2))
  expect_that(top$itemCount > 4 * 1000 * 1000, is_true())
  
  # Save the number of media in all iDigBio for later tests
  #all_count <- top$itemCount
})

test_that("specifying field and number of tops works", {
  testthat::skip_on_cran()
  top <- idig_top_media(top_fields=c(field), count=count)
  
  expect_that(top, is_a("list"))
  expect_that(length(top[[field]]), equals(count))
  expect_that(top[[field]][[most]][["itemCount"]] > 10 * 1000, is_true())
  
  # Deprecating this since Alex changed the erorr behavior to tell you when
  # JSON is bad or field unknown, no longer just spits out all iDigBio
  # Still looking at all of iDigBio, assume things are not changing too fast
  #expect_that(abs(top$itemCount - all_count) < 1000, is_true())
})

test_that("searches work", {
  testthat::skip_on_cran()
  top <- idig_top_media(mq=list("dqs"=dqs), top_fields=c(field),
                         count=count)
  
  expect_that(top, is_a("list"))
  #expect_that(top$itemCount < all_count, is_true())
  expect_that(top$itemCount > 0, is_true())
})


# WIP

# Save the number of genus records for later tests
#genus_count <- top$itemCount

# Multiple fields
#top <- idig_top_media(rq=list("genus"=genus), top_fields=c(field, 
#                       "scientificname"), count=count)
#expect_that(top, is_a("list"))
#expect_that(abs(top$itemCount - genus_count) < 100, is_true())
#expect_that(top[[field]][[most]][["scientificname"]][[scientificname]][["itemCount"]]
#            > 1000, is_true())