context("Reshaping data")

test_that("fast reshape I", {    
    m <- lvm()
    regression(m,c(y1,y2,y3)~x) <- c(1,10,100)
    distribution(m,~x) <- f <- function(n,...) rbinom(n,1,0.5)+1
    d <- sim(m,10); 
    dd <- fast.reshape(d,var="y")
    d1 <- fast.reshape(dd,id="id")
    expect_true(sum((d[,endogenous(m)]-d1[,endogenous(m)])^2)<1e-20)
    d2 <- fast.reshape(dd,id="id",var="y",num="num")
    expect_true(sum((d-d2[,colnames(d)])^2)<1e-20)
})


test_that("fast reshape II", {
    testdata <- data.frame(hour=c(12,13,14,11,12,14,15,16),id=c(1,1,1,2,2,3,3,3),y=round(rnorm(8),2))
    widetest <- reshape(testdata,v.names="y",idvar="id",direction="wide",timevar="hour")
    wide <- fast.reshape(testdata,varying="y",id="id",num="hour",sep=".")
    expect_equivalent(widetest,wide[,colnames(widetest)])
})


## fast.reshape(fast.reshape(d,var=c("y","z","w")),id="id",var=c("y","z","w"))
## library(mets)
## x <- matrix(1:10,5,2)
## x[3,2] <- 8
## x[3,2] <- NA
## cluster <- c(1,1,2,2,3)
## x <- cbind(x,cluster)
## x
## ud <- fast.reshape(data.frame(x),"cluster")
## ud
## ###
## out=cluster.index(cluster)
## out
## ###
## ud <- faster.reshape(x,cluster)
## ud
## ud <- faster.reshape(data.frame(x),cluster)
## ud
## ###
## colnames(x) <- c("y1","y2","cluster")
## x
## ud <- fast.reshape(data.frame(x),"cluster")
## ud
## ###
## num <- c(2,1,1,2,2)
## x
## out <- faster.reshape(x,cluster)
## out
## out <- faster.reshape(x,cluster,num=num)
## out

