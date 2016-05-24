context("Test for lets_transf")

status <- sample(c("EN","VU", "NT", "CR", "DD", "LC"), 30, replace=TRUE)
TE <- "Threatened"
NT <- "Non-Threatened"
new2 <- c(TE, TE, NT, TE, "Data Deficient", NT)
old <- c("EN","VU", "NT", "CR", "DD", "LC")
new3 <- c(1, 2, 1, 2, 3, 40)

test_that("lets_transf works fine", {
  
  statustrans <- lets.transf(status, old, new2, NUMERIC = FALSE)
  expect_true(is.character(statustrans))
    
  statustrans2 <- lets.transf(status, old, new3, NUMERIC = TRUE)
  expect_true(is.numeric(statustrans2))
  
})
