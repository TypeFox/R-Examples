context("get.theta")
test_that("get.theta returns Inf when all relations are 1", { #Would throwing an error be better?
    #generate a set of 100 random points even labeled between the two
    x<-cbind(rep(c(1,2),50), x=runif(100,0,100), y=runif(100,0,100))

    colnames(x) <-c("type","x","y")

    test <- function(a,b) {return(1)}

    #with no lower limit
    res <- get.theta(x,test,seq(10,100,10))
    expect_that(res,equals(rep(Inf,10)))

    #with lower and upper limit
    res <- get.theta(x,test,seq(10,100,10), seq(0,90,10))
    expect_that(res,equals(rep(Inf,10)))
})


test_that("get.theta returns 0 when all relations are 2", {
     #generate a set of 100 random points even labeled between the two
    x<-cbind(rep(c(1,2),50), x=runif(100,0,100), y=runif(100,0,100))

    colnames(x) <-c("type","x","y")

    test <- function(a,b) {return(2)}

    #with no lower limit
    res <- get.theta(x,test,seq(10,100,10))
    expect_that(res,equals(rep(0,10)))

    #with lower and upper limit
    res <- get.theta(x,test,seq(10,100,10), seq(0,90,10))
    expect_that(res,equals(rep(0,10)))
})

test_that("get.theta returns appropriate values for cannonical test case 1 (equilateral triangle)", {

    #in equilateral triangle case, exactly half of the pairs
    #or interest are 1,1...meaning that the odds should be 1 in
    #each case

    x <- rbind(c(1,0,0), c(1,1,0),c(2,.5,sqrt(.75)))
    colnames(x) <-c("type","x","y")

    test <- function(a,b) {
        if (a[1] != 1) return(3)
        if (b[1] == 2) return(1)
        return(2)
    }


    #first no lower limit
    res <- get.theta(x,test,1.5)
    res2 <- get.theta.typed(x,1,2,1.5)

    expect_that(res, equals(1))
    expect_that(res2, equals(1))

    #now with a lower limit

    res <- get.theta(x,test,1.5,.5)
    res2 <- get.theta.typed(x,1,2,1.5,.5)

    expect_that(res, equals(1))
    expect_that(res2, equals(1))

})

test_that("get.theta returns appropriate values cannonical test case 2 (points on a line)", {
    x<-rbind(c(1,0,0), c(2,1,0), c(2,-1,0), c(3,2,0),
             c(2,-2,0), c(3,3,0),c(3,-3,0))

    colnames(x) <-c("type","x","y")

    test <- function(a,b) {
        if (a[1] != 1) return(3)
        if (b[1] == 2) return(1)
        return(2)
    }


    #pi 0,1.5 should be 1, 1.5-2.5 should be 0.5 and 2.5+ should be 0
    res <- get.theta(x, test, c(1.5,2.5,Inf), c(0,1.5,2.5))
    res2 <- get.theta.typed(x, 1, 2, c(1.5,2.5,1000), c(0,1.5,2.5))

    expect_that(res,equals(c(Inf,1,0)))
    expect_that(res2,equals(c(Inf,1,0)))

})

test_that("get.theta and get.theta.typed have same behavior on random data", {

    #generate a set of 1000 random points even labeled between the two
    x<-cbind(rep(c(1,2),50), x=runif(100,0,100), y=runif(100,0,100))

    colnames(x) <-c("type","x","y")

    test <- function(a,b) {
        if (a[1] != 1) return(3)
        if (b[1] == 1) return(1)
        return(2)
    }

    #no lower limit
    res1 <- get.theta(x,test,seq(10,100,10))
    res2 <- get.theta.typed(x, 1,1, seq(10,100,10))
    expect_that(res1,equals(res2))

    #lower limit
    res1 <- get.theta(x,test,seq(10,100,10), seq(0,90,10))
    res2 <- get.theta.typed(x, 1,1, seq(10,100,10), seq(0,90,10))
    expect_that(res1,equals(res2))
})

test_that("get.theta returns identical results regardless of column order",
          {
              x<-cbind(rep(c(1,2),50), x=runif(100,0,100),
                       y=runif(100,0,100))

              colnames(x) <-c("type","x","y")

              test <- function(a,b) {
                  if (a[1] != 1) return(3)
                  if (b[1] == 1) return(1)
                  return(2)
              }

              res1 <- get.theta(x,test,seq(10,100,10), seq(0,90,10))

              test <- function(a,b) {
                  if (a[3] != 1) return(3)
                  if (b[3] == 1) return(1)
                  return(2)
              }

              res2 <- get.theta(x[,c(3,2,1)],test,seq(10,100,10), seq(0,90,10))

              test <- function(a,b) {
                  if (a[2] != 1) return(3)
                  if (b[2] == 1) return(1)
                  return(2)
              }

              res3 <- get.theta(x[,c(2,1,3)],test,seq(10,100,10), seq(0,90,10))

              expect_that(res1, equals(res2))
              expect_that(res2, equals(res3))

          })

test_that ("get.theta fails nicely if x and y column names are not provided", {
    x<-cbind(rep(c(1,2),500), a=runif(1000,0,100), b=runif(1000,0,100))

    test <- function(a,b) {
        if (a[1] != 2) return(3)
        if (b[1] == 3) return(1)
        return(2)
    }

    expect_that(get.theta(x,test,seq(10,50,10), seq(0,40,10)),
                throws_error("unique x and y columns must be defined"))

})
