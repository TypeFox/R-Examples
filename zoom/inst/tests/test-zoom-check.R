test_that("multipancPoint calculates fine new limits",{
expect_equal(multipancPoint(c(0,1),fact=2,NULL,point=0.1),c(0.05,0.55))
expect_equal(multipancPoint(c(-1,1),fact=1,NULL,point=0),c(-1,1))
expect_equal(multipancPoint(c(-1,1),fact=0.5,NULL,point=0),c(-2,2))
})


