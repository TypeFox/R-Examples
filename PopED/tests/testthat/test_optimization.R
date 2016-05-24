context("Optimization")

test_that("optim_ARS() works", {
  
  ex_string <- ex_to_string("examples_fcn_doc/examples_optim_ARS.R",comment_dontrun=comment_dontrun)
  sink("tmp.txt")
  eval(parse(text=ex_string))
  sink()
  file.remove("tmp.txt")
  
  expect_true(res1$ofv <= 159.0012) # check that things are minimizing
  
  # check that box constraints hold
  expect_true(all(res_box$par<=4)) 
  expect_true(all(res_box$par >= 2))
              
  # check that there are only these allowed values in the solution
  expect_true(res_int$par %in% seq(-50,100,by=1))
  
  # check that function in maximizing
  expect_true(res_max$ofv >= 0)
  
})