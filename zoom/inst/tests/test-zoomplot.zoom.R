library("testthat")
## need to comment the following for the non interactive device
## based tests of the CRAN
test_that("test basic zoomplot.zoom fact is ok",{
plot(rnorm(100),rnorm(100),xlim=c(-2,2),ylim=c(-2,2))
zoomplot.zoom(fact=2)
newLim<-par("usr")
absDiffx<-abs(diff(newLim[1:2]))
# cat("newLim:",newLim,"absdiffx:",absDiffx,"\n")
expect_true(absDiffx<3)
expect_true(absDiffx>2)
dev.off()
})

test_that("test zoomplot.zoom direct xlim, ylim set is ok",{
plot(rnorm(100),rnorm(100),xlim=c(-1,1),ylim=c(-1,1))
zoomplot.zoom(xlim=c(-2,2),ylim=c(-2,2))
newLim<-par("usr")
absDiffx<-abs(diff(newLim[1:2]))
# cat("newLim:",newLim,"absdiffx:",absDiffx,"\n")
expect_true(absDiffx>4)
expect_true(absDiffx<5)
dev.off()
})

test_that("test zoomplot.zoom fact + stable pointer works",{
plot(rnorm(100),rnorm(100),xlim=c(-2,2),ylim=c(-2,2))
zoomplot.zoom(fact=2,x=0.5,y=0.5)
newLim<-par("usr")
absDiffx<-abs(diff(newLim[1:2]))
absDiffy<-abs(diff(newLim[3:4]))
# cat("newLim:",newLim,"absdiffx:",absDiffx,"\n")
expect_true(absDiffx<2.5)
expect_true(absDiffx>2)
expect_equal(absDiffx,absDiffy)
expect_true(abs(grconvertX(0.5,"user","ndc")-0.6)<0.1)
expect_true(abs(grconvertY(0.5,"user","ndc")-0.6)<0.1)
dev.off()
})

