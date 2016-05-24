context("Test for iucn_cont")

status <- sample(c("EN","VU", "NT", "CR", "DD", "LC", "EX", "NE"),
                 30, replace = TRUE)
data(IUCN)


test_that("iucn_cont works fine", {
  
  transV <- lets.iucncont(status)
  expect_true(is.numeric(transV))
  
  #matrix transformation
  transM <- lets.iucncont(IUCN)  
  expect_true(is.data.frame(transM))
  
  transV2 <- lets.iucncont(status, dd = 2, ne = 1)
  expect_true(is.numeric(transV2))
  
  
})
