# library(testthat)
# library(ff)
# 
# test_that("risplit works",{
#    x <- 1:10
#    f <- as.factor(rep(c("M","F"), 5))
#    risplit(x, f)
# })
# 
# test_that("risplit works for ff",{
#    x <- 1:10
#    fa <- as.factor(rep(c("M","F"), 5))
#    
#    xf <- ff(x)
#    faf <- ff(fa)
#    
#    rs1 <- risplit(x,fa)
#    rs2 <- risplit(xf,faf)
#    
#    expect_equivalent( rs1
#                     , rs2
#                     )
# 
# })
# 
# test_that("risplit works for ffdf",{
#    xdf <- as.ffdf(airquality)
#    g <- xdf$Month
#    #l <- risplit(xdf, g)
# })