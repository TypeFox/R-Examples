context("GPCM")


THRES  <- matrix(c(-2,-1.23,1.11,3.48,1
                   ,2,-1,-0.2,0.5,1.3,-0.8,1.5),nrow=2)

sl     <- c(0.5,1,1.5,1.1,1,0.98)
awmatrix <- matrix(c(1,0,2,0,1,1,1,0,0,1
                     ,2,0,0,0,0,0,0,0,0,1,1,2,2,1,1,1,1,0,0,1),byrow=TRUE,nrow=5)


awmatrix_2 <- awmatrix
awmatrix_2[2,1] <- 3

THRESna  <- THRES

THRESna[,1] <- NA

test_that("errors and warnings are thrown in the gpcm when necessary",{
  expect_that(PP_gpcm(respm = awmatrix,thres = THRES, slopes = sl,type = "robust"),gives_warning())
  expect_that(PP_gpcm(respm = awmatrix,thres = THRES, slopes = sl,type = "A"),throws_error())
  expect_that(PP_gpcm(respm = awmatrix[,-1],thres = THRES, slopes = sl),throws_error())
  expect_that(PP_gpcm(respm = awmatrix[,-1],thres = THRES[,-1], slopes = sl),throws_error())
  expect_that(PP_gpcm(respm = awmatrix[,-1],thres = THRES, slopes = sl[-1]),throws_error())
  expect_that(PP_gpcm(respm = awmatrix,thres = THRESna, slopes = sl),throws_error())
  expect_that(PP_gpcm(respm = awmatrix_2,thres = THRES, slopes = sl),throws_error())
  expect_that(PP_gpcm(respm = awmatrix,thres = THRES, slopes = sl,theta_start = 2),throws_error())
})

