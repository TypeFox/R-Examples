test_that("Read espresso file", {
   data(l.small)
   esp1 <- logicopt(l.small,4,3)
   file <- system.file("extdata/espresso/small.esp", package="LogicOpt")
   esp7 <- logicopt(esp_file=file)
   expect_equal(esp1,esp7)
  
   expect_warning(logicopt(esp_file=file,find_dc=TRUE),
"When esp_file is used, options find_dc, n_in, n_out, and input_sizes are ignored.")

   expect_warning(logicopt(esp_file=file,n_in=4),
"When esp_file is used, options find_dc, n_in, n_out, and input_sizes are ignored.")

   expect_warning(logicopt(esp_file=file,n_out=3),
"When esp_file is used, options find_dc, n_in, n_out, and input_sizes are ignored.")

   expect_warning(logicopt(esp_file=file,input_sizes=c(2,2,2,2)),
"When esp_file is used, options find_dc, n_in, n_out, and input_sizes are ignored.")

})

test_that("Logicopt echo tests", {
   file <- system.file("extdata/espresso/robot1_in.esp", package="LogicOpt")
   robot1_echo <- logicopt(esp_file=file,mode="echo")
   robot1_esp <- logicopt(esp_file=file,mode="espresso")
   robot1_opt <- logicopt(robot1_echo[[1]],8,3,mode="espresso")
   expect_equal(robot1_esp,robot1_opt)

   file <- system.file("extdata/espresso/small2.esp", package="LogicOpt")
   small <- logicopt(esp_file=file,mode="echo")
   small_pr1 <- logicopt(small[[1]],4,3,mode="primes")
   small_pr2 <- logicopt(esp_file=file,mode="primes")
   expect_equal(small_pr1,small_pr2)

})
