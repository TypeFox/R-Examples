context("Clustered survival data")
test_that("clustersruv",{
    library(prodlim)
    ## if (!is.function("cluster")) cluster <- function(x)x
    clusterTestData <- data.frame(midtimeX=1:8,eventX=c(0,"pn","pn",0,0,0,0,0),patientid=c(1,1,2,2,3,3,4,4),AnyCrownFracture=c(1,1,1,1,2,2,2,2))
    a <- prodlim(Hist(midtimeX,eventX=="pn")~cluster(patientid)+AnyCrownFracture,data=clusterTestData)
    b <- prodlim(Hist(midtimeX,eventX=="pn")~cluster(patientid),data=clusterTestData[clusterTestData$AnyCrownFracture==1,])
    c <- prodlim(Hist(midtimeX,eventX=="pn")~cluster(patientid),data=clusterTestData,subset=clusterTestData$AnyCrownFracture==1)
    d <- prodlim(Hist(midtimeX,eventX=="pn")~1,data=clusterTestData[clusterTestData$AnyCrownFracture==2,])
    expect_equal(round(as.numeric(summary(a)$table[[1]][,c("se.surv")]),5),c(0,0.20951,0.10476,0.10476,NA,NA,NA,NA))
    expect_equal(summary(b), summary(c))
})
