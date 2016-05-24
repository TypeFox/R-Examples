test_that("espresso errors", {
   data(l.small)
   esp1 <- logicopt(l.small,4,3)
  
   expect_error(logicopt(esp1[[1]],4,4),"is not equal to")
   expect_error(logicopt(esp1[[1]],5,3),"is not equal to")
   expect_error(logicopt(esp_file="somerandomfile"),"does not exist")
   expect_error(logicopt(esp1[[1]],4,3,TRUE),"espresso failed")
 
   file <- system.file("extdata/espresso/error.esp", package="LogicOpt")
   #espresso messaging messes up test()
   #expect_error(logicopt(esp_file=file),"espresso failed")

})

