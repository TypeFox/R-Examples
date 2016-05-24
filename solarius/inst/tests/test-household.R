
context("household")

test_that("house-hold default (NA)", {
  data(dat30)
  
  dat30 <- mutate(dat30,
    hhid = famid)  
    
  mod <- solarPolygenic(trait1 ~ 1, dat30)

  expect_true(any(grepl("C2 parameter", mod$solar$files$model$polygenic.out)))
})

test_that("house-hold TRUE", {
  data(dat30)
  
  dat30 <- mutate(dat30,
    hhid = famid)  
    
  mod <- solarPolygenic(trait1 ~ 1, dat30, household = TRUE)

  expect_true("c2" %in% mod$vcf$varcomp)
})

