
context("4PL model")

# ------------------------- testing 1>>>


THRESx  <- c(-2,-1,1,2)
sl     <- c(0.5,1,1.5,1.1)

awmatrix <- matrix(c(1,0,1,0,1,1,1,0),nrow=2,byrow=TRUE)
awmatrix <- rbind(awmatrix,c(1,1,1,1),c(0,0,0,0),c(1,0,1,0),c(1,0,1,0),c(1,1,0,0))

ua     <- c(0.98,0.85,0.9,0.95)
la     <- c(0,0.05,0.12,0.001)


estmod <- rep(c("mle","wle","map"),4)
LA <- vector(mode="list",length=length(estmod))
UA <- vector(mode="list",length=length(estmod))
for(i in c(4,5,6,10,11,12))
  {
    LA[[i]] <- la      
  }

for(i in c(7,8,9,10,11,12))
  {
    UA[[i]] <- ua
  }
    

res234pl_dup1 <- vector(mode="list",length=length(estmod))
res234pl_dup2 <- vector(mode="list",length=length(estmod))


for(a in 1:length(estmod))
  {
    res234pl_dup1[[a]]  <- PP_4pl(awmatrix,THRESx,slopes = sl,type = estmod[[a]],ctrl = list(killdupli=TRUE),upperA = UA[[a]],lowerA = LA[[a]])
  
    res234pl_dup2[[a]] <- PP_4pl(awmatrix,THRESx,slopes = sl,type = estmod[[a]],ctrl = list(killdupli=FALSE),upperA = UA[[a]],lowerA = LA[[a]])
    
  }


#t
test_that("Output = the same - with or without removing duplicates",{

for(te in 1:length(estmod))
{
expect_that(res234pl_dup1[[te]],equals(res234pl_dup1[[te]]))  
}

})




# ------------------------- testing 2>>>


set.seed(1523)
# intercepts
diffpar <- seq(-3,3,length=12)
la     <- round(runif(12,0,0.25),2)
ua     <- round(runif(12,0.8,1),2)
# slope parameters
sl     <- round(runif(12,0.5,1.5),2)

# antwortmatrix (für neue MLE routine)
awm <- matrix(sample(0:1,10*12,replace=TRUE),ncol=12)
awm2 <- awm
awm2[3,2] <- 2

diffparM <- rbind(0,diffpar)
diffparM2 <- diffparM
diffparM2[2,4] <- NA

b1 <- PP_4pl(respm = awm,thres = diffpar, slopes = sl,type = "wle")
b2 <- PP_4pl(respm = awm,thres = diffparM, slopes = sl,type = "wle")

b3 <- PP_4pl(respm = awm,thres = diffpar, slopes = sl,type = "mle")
b4 <- PP_4pl(respm = awm,thres = diffparM, slopes = sl,type = "mle")

#t
test_that("Output = the same on 1,2,3,4pl - with vector or matrix input",{
  expect_that(b1[[1]],equals(b2[[1]]))
  expect_that(b3[[1]],equals(b4[[1]]))
})



# ------------------------- testing 3>>>

#t
test_that("errors - warnings misspelling and length #1",{
  expect_that(PP_4pl(respm = awm,thres = diffpar, slopes = sl,type = "aaa"), throws_error())
  expect_that(PP_4pl(respm = awm,thres = diffpar[-1], slopes = sl), throws_error())
  expect_that(PP_4pl(respm = awm[,-1],thres = diffpar, slopes = sl), throws_error()) 
  expect_that(PP_4pl(respm = awm,thres = diffpar[-1], slopes = sl[-1]), throws_error()) 
  expect_that(PP_4pl(respm = awm,thres = diffpar, slopes = sl,lowerA = la[-1]), throws_error()) 
  expect_that(PP_4pl(respm = awm,thres = diffpar[-1], slopes = sl[-1],upperA = ua[-1]), throws_error()) 
  expect_that(PP_4pl(respm = awm,thres = diffpar[-1], slopes = sl[-1],upperA = ua[-1],lowerA = la[-1]), throws_error())
  expect_that(PP_4pl(respm = awm,thres = diffpar[-1], slopes = sl[-1],upperA = ua[-1],lowerA = la[-1]), throws_error())
  expect_that(PP_4pl(respm = awm2,thres = diffpar, slopes = sl), throws_error())
  expect_that(PP_4pl(respm = awm2,thres = diffparM2, slopes = sl), throws_error())
})



















