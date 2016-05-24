context("Concatenating")

test_that("concatenate editmatrix",{
    
    E1  <- editmatrix("x1 + x2 > 0")
    E2  <- editmatrix("x2 + x3 < 10")
    
    E <- editmatrix(c("x1 + x2 > 0", "x2 + x3 < 10"))
    #print(c.editmatrix(E1,E2))
    #print(E)
    expect_equivalent(E, c(E1,E2))
    e <- editmatrix(expression())
    expect_equivalent(nedits(c(e,e)),0)
    expect_equivalent(nedits(c(e,NULL)),0)
})

test_that("concatenate editmatrix",{  
  E <- editmatrix(c("x1 + x2 > 0", "x2 + x3 < 10"))
  expect_equivalent(E, c(E))
})


