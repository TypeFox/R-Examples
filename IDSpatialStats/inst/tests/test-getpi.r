context("get.pi")
test_that("get.pi returns 1 when labels are ignored", {

    #generate a set of 100 random points even labeled between the two
    x<-cbind(rep(c(1,2),50), x=runif(100,0,100), y=runif(100,0,100))

    colnames(x) <-c("type","x","y")

    test <- function(a,b) {return(1)}

    #with no lower limit
    res <- get.pi(x,test,seq(10,100,10))
    expect_that(res,equals(rep(1,10)))

    #with lower and upper limit
    res <- get.pi(x,test,seq(10,100,10), seq(0,90,10))
    expect_that(res,equals(rep(1,10)))
})

test_that("get.pi returns appropriate values for cannonical test case 1 (equilateral triangle)", {

    x <- rbind(c(1,0,0), c(1,1,0),c(2,.5,sqrt(.75)))
    colnames(x) <-c("type","x","y")

    test <- function(a,b) {
        if (a[1] != 1) return(3)
        if (b[1] == 2) return(1)
        return(2)
    }

    #first no lower limit
    res <- get.pi(x,test,1.5)
    res2 <- get.pi.typed(x,1,2,1.5)

    expect_that(res, equals(.5))
    expect_that(res2, equals(.5))

    #now with a lower limit

    res <- get.pi(x,test,1.5,.5)
    res2 <- get.pi.typed(x,1,2,1.5,.5)

    expect_that(res, equals(.5))
    expect_that(res2, equals(.5))

})

test_that("get.pi returns appropriate values cannonical test case 2 (points on a line)", {
    x<-rbind(c(1,0,0), c(2,1,0), c(2,-1,0), c(3,2,0),
             c(2,-2,0), c(3,3,0),c(3,-3,0))

    colnames(x) <-c("type","x","y")

    test <- function(a,b) {
        if (a[1] != 1) return(3)
        if (b[1] == 2) return(1)
        return(2)
    }

    #pi 0,1.5 should be 1, 1.5-2.5 should be 0.5 and 2.5+ should be 0
    res <- get.pi(x, test, c(1.5,2.5,Inf), c(0,1.5,2.5))
    res2 <- get.pi.typed(x, 1, 2, c(1.5,2.5,1000), c(0,1.5,2.5))

    expect_that(res,equals(c(1,0.5,0)))
    expect_that(res2,equals(c(1,0.5,0)))

})


test_that("get.pi and get.pi.typed have same behavior on random data", {

    #generate a set of 1000 random points even labeled between the two
    x<-cbind(rep(c(1,2),50), x=runif(100,0,100), y=runif(100,0,100))

    colnames(x) <-c("type","x","y")

    test <- function(a,b) {
        if (a[1] != 1) return(3)
        if (b[1] == 1) return(1)
        return(2)
    }

    #no lower limit
    res1 <- get.pi(x,test,seq(10,100,10))
    res2 <- get.pi.typed(x, 1,1, seq(10,100,10))
    expect_that(res1,equals(res2))

    #lower limit
    res1 <- get.pi(x,test,seq(10,100,10), seq(0,90,10))
    res2 <- get.pi.typed(x, 1,1, seq(10,100,10), seq(0,90,10))
    expect_that(res1,equals(res2))
})

test_that("get.pi returns identical results regardless of column order",
          {
              x<-cbind(rep(c(1,2),50), x=runif(100,0,100),
                       y=runif(100,0,100))

              colnames(x) <-c("type","x","y")

              test <- function(a,b) {
                  if (a[1] != 1) return(3)
                  if (b[1] == 1) return(1)
                  return(2)
              }

              res1 <- get.pi(x,test,seq(10,100,10), seq(0,90,10))

              test <- function(a,b) {
                  if (a[3] != 1) return(3)
                  if (b[3] == 1) return(1)
                  return(2)
              }

              res2 <- get.pi(x[,c(3,2,1)],test,seq(10,100,10), seq(0,90,10))

              test <- function(a,b) {
                  if (a[2] != 1) return(3)
                  if (b[2] == 1) return(1)
                  return(2)
              }

              res3 <- get.pi(x[,c(2,1,3)],test,seq(10,100,10), seq(0,90,10))

              expect_that(res1, equals(res2))
              expect_that(res2, equals(res3))

          })


test_that ("get.pi fails nicely if x and y column names are not provided", {
    x<-cbind(rep(c(1,2),500), a=runif(1000,0,100), b=runif(1000,0,100))

    test <- function(a,b) {
        if (a[1] != 2) return(3)
        if (b[1] == 3) return(1)
        return(2)
    }

    expect_that(get.pi(x,test,seq(10,50,10), seq(0,40,10)),
                throws_error("unique x and y columns must be defined"))

})


##################DEPRECATED TESTS...TAKE TO LONG...NOW USING SMALLER CANONICAL
##################TESTS THAT HAVE VALUES THAT CAN BE WORKED OUT BY HAND
## test_that("get.pi and get.pi.typed have same behavior and both return about .5 when they should", {

##     set.seed(787)

##     #generate a set of 1000 random points even labeled between the two
##     x<-cbind(rep(c(1,2),500), x=runif(1000,0,100), y=runif(1000,0,100))

##     colnames(x) <-c("type","x","y")

##     test <- function(a,b) {
##         if (a[1] != 1) return(3)
##         if (b[1] == 1) return(1)
##         return(2)
##     }

##     #no lower limit
##     res1 <- get.pi(x,test,seq(10,100,10))
##     res2 <- get.pi.typed(x, 1,1, seq(10,100,10))
##     expect_that(res1,equals(res2))
##     expect_that(round(res1,1),equals(rep(.5,10)))
##     expect_that(round(res2,1),equals(rep(.5,10)))


##     #lower limit
##     res1 <- get.pi(x,test,seq(10,100,10), seq(0,90,10))
##     res2 <- get.pi.typed(x, 1,1, seq(10,100,10), seq(0,90,10))
##     expect_that(res1,equals(res2))
##     expect_that(round(res1,1),equals(rep(.5,10)))
##     expect_that(round(res2,1),equals(rep(.5,10)))
## })


## test_that("get.pi with equality test returns about .5 when points are uniform",
##       {
##           set.seed(787)

##           #generate a set of 1000 random points even labeled between the two
##           x<-cbind(rep(c(1,2),500), x=runif(1000,0,100), y=runif(1000,0,100))

##           colnames(x) <-c("type","x","y")

##           test <- function(a,b) {
##               if (a[1]==b[1]) return(1)
##               return(2)
##           }

##           #no lower limit
##           res1 <- get.pi(x,test,seq(10,100,10))
##           expect_that(round(res1,1),equals(rep(.5,10)))

##           #lower limit
##           res1 <- get.pi(x,test,seq(10,100,10), seq(0,90,10))
##           expect_that(round(res1,1),equals(rep(.5,10)))
##       })


## test_that("get.pi is montonically decreasing for normally distributed clusters",
##           {
##               set.seed(787)

##               #first generate 500 random uniform points
##               x<-cbind(1, x=runif(500,0,100), y=runif(500,0,100))
##               colnames(x) <-c("type","x","y")

##               #add a seed point
##               x<-rbind(x,c(2,50,50))

##               #generate 500 normally distibuted points around this
##               x<-rbind(x,cbind(3,rnorm(1000,50,20),rnorm(1000,50,20)))

##               #check wit get.pi.typed
##               res1 <- get.pi.typed(x,2,3,seq(10,50,10), seq(0,40,10))

##               expect_that(res1[1]>res1[2] & res1[2]>res1[3] &
##                           res1[3]>res1[4] & res1[4]>res1[5], is_true())


##               #do test for not pi version
##               test <- function(a,b) {
##                   if (a[1] != 2) return(3)
##                   if (b[1] == 3) return(1)
##                   return(2)
##               }

##               res2 <- get.pi(x,test,seq(10,50,10), seq(0,40,10))


##               expect_that(res2[1]>res2[2] & res2[2]>res2[3] &
##                           res2[3]>res2[4] & res2[4]>res2[5], is_true())
##               expect_that(res1,equals(res2))
##           })

## test_that("get.pi returns identical results regardless of column order",
##           {
##               set.seed(787)

##               x<-cbind(rep(c(1,2),500), x=runif(1000,0,100), y=runif(1000,0,100))
##               colnames(x) <-c("type","x","y")

##               test <- function(a,b) {
##                   if (a[1] != 1) return(3)
##                   if (b[1] == 1) return(1)
##                   return(2)
##               }

##               res1 <- get.pi(x,test,seq(10,100,10), seq(0,90,10))

##               test <- function(a,b) {
##                   if (a[3] != 1) return(3)
##                   if (b[3] == 1) return(1)
##                   return(2)
##               }

##               res2 <- get.pi(x[,c(3,2,1)],test,seq(10,100,10), seq(0,90,10))

##               test <- function(a,b) {
##                   if (a[2] != 1) return(3)
##                   if (b[2] == 1) return(1)
##                   return(2)
##               }

##               res3 <- get.pi(x[,c(2,1,3)],test,seq(10,100,10), seq(0,90,10))

##               expect_that(res1, equals(res2))
##               expect_that(res2, equals(res3))

##           })



