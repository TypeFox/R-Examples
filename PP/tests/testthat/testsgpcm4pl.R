
context("GPCM&4PL model")

# some threshold parameters
THRES  <- matrix(c(-2,-1.23,1.11,3.48,1
                   ,2,-1,-0.2,0.5,1.3,-0.8,1.5),nrow=2)
# slopes
sl     <- c(0.5,1,1.5,1.1,1,0.98)

THRESx <- THRES
THRESx[2,1:3] <- NA
THRESx2 <- THRESx

THRESx2[,4] <- NA

#
la     <- c(0.02,0.1,0,0,0,0)
ua     <- c(0.97,0.91,1,1,1,1)

awm <- matrix(c(1,0,1,0,1,1,1,0,0,1
                     ,2,0,0,0,0,0,0,0,0,1
                     ,0,2,2,1,1,1,1,0,0,1),byrow=TRUE,nrow=5)

# create model2est
# this function tries to help finding the appropriate 
# model by inspecting the THRESx.
model2est <- findmodel(THRESx)

awm3 <- awm2 <- awm 
awm2[3,1] <- 3

#awm3[,1] <- NA

# ------------------------- testing 1>>>

#t
test_that("errors - warnings misspelling and length #1",{
  expect_that(PPall(respm = awm,thres = THRESx, slopes = sl,lowerA=la,upperA=ua,type = "aaa",model2est=model2est), throws_error())
  
  expect_that(PPall(respm = awm,thres = THRESx[,-1], slopes = sl,lowerA=la,upperA=ua,model2est=model2est), throws_error())
  
  expect_that(PPall(respm = awm,thres = THRESx, slopes = sl[-1],lowerA=la,upperA=ua,model2est=model2est), throws_error())
  
  expect_that(PPall(respm = awm,thres = THRESx, slopes = sl,lowerA=la[-1],upperA=ua,model2est=model2est), throws_error())
  
  expect_that(PPall(respm = awm,thres = THRESx, slopes = sl,lowerA=la,upperA=ua[-1],model2est=model2est), throws_error())
  
  expect_that(PPall(respm = awm[,-1],thres = THRESx, slopes = sl,lowerA=la,upperA=ua,model2est=model2est), throws_error())
  
  expect_that(PPall(respm = awm2,thres = THRESx, slopes = sl,lowerA=la,upperA=ua,model2est=model2est), throws_error())
  
  expect_that(PPall(respm = awm2,thres = THRESx2, slopes = sl,lowerA=la,upperA=ua,model2est=model2est), throws_error())
  
  expect_that(PPall(respm = awm,thres = THRESx, slopes = sl,lowerA=la,upperA=ua,model2est=model2est, type="robust"), gives_warning())
})
