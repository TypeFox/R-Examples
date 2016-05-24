#R

context("Test 8 reps data")

test_that("get angle of two map regions", {
  
  
  
	expect_true( recmap:::get_angle(x0=0, y0=0, x1=0, y1=1)  / pi == 0)
	expect_true(recmap:::get_angle(x0=0, y0=0, x1=0, y1=-1)  / pi == 1 )
	expect_true(recmap:::get_angle(x0=0, y0=0, x1=1, y1=-1)  / pi == 0.75 )
       expect_true( recmap:::get_angle(x0=0, y0=0, x1=-1, y1=-1)  / pi == -0.75 )
       expect_true( recmap:::get_angle(x0=0, y0=0, x1=-1, y1=1)  / pi == -0.25 )
})
