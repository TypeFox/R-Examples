context("Summary")
test_that("Summary", {
  suppressWarnings({
    test <- defineClass("test", {
      x <- publicValue(1)
      .y <- NULL
      doSomething <- publicFunction(function() {
        .y <<- .y + 1
        invisible(self)
      })
    })
  })
  
  instance <- test()
  df <- summary(instance)
  expect_is(df, "data.frame")
  expect_equal(nrow(df), 7)
  expect_equal(sort(df$Name), 
               sort(c("doSomething", "self", ".self", "x", "x.validity", "x.x", ".y")))
  
})