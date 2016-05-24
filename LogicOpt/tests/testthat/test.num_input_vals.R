test_that("Find unique input values", {
   data(l.partybans.1)
   vals <- num_input_values(l.partybans.1,5)
   expect_equal(vals,c(3,3,3,3,2))

   data(l.small)
   vals <- num_input_values(l.small,4)
   expect_equal(vals,c(2,2,2,2))
})
