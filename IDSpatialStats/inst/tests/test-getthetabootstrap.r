context("get.theta.bootstrap")

test_that("get.theta.boostrap runs and returns Inf when all relations are 1", {

    x<-cbind(rep(c(1,2),50), x=runif(100,0,100), y=runif(100,0,100))

    colnames(x) <-c("type","x","y")

    test <- function(a,b) {return(1)}

    #should return a matrix of all ones
    res <- get.theta.bootstrap(x, test, seq(10,100,10), seq(0,90,10), 20)
    expect_that(sum(!is.infinite(res)),equals(0))
    expect_that(nrow(res),equals(20))

})


test_that("get.theta.ci returns bootstrap cis when same seed", {
    x<-cbind(rep(c(1,2),50), x=runif(100,0,100), y=runif(100,0,100))

    colnames(x) <-c("type","x","y")

    test <- function(a,b) {
        if (a[1] != 1) return(3)
        if (b[1] == 2) return(1)
        return(2)
    }

    set.seed(787)
    res <- get.theta.bootstrap(x, test, seq(15,45,15), seq(0,30,15), 20)

    set.seed(787)
    ci1 <- get.theta.ci(x, test, seq(15,45,15), seq(0,30,15), 20)

    expect_that(as.numeric(ci1[,1]),
                equals(as.numeric(quantile(res[,1],
                                           probs=c(.025,.975),
                                           na.rm=T))))

    expect_that(as.numeric(ci1[,2]),
                equals(as.numeric(quantile(res[,2],
                                           probs=c(.025,.975),
                                           na.rm=T))))

    expect_that(as.numeric(ci1[,3]),
                equals(as.numeric(quantile(res[,3],
                                           probs=c(.025,.975),
                                           na.rm=T))))

})


test_that("performs correctly for test case 1 (equilateral triangle)", {
    x <- rbind(c(1,0,0), c(1,1,0),c(2,.5,sqrt(.75)))
    colnames(x) <-c("type","x","y")

    test <- function(a,b) {
        if (a[1] != 1) return(3)
        if (b[1] == 2) return(1)
        return(2)
    }

    res <- get.theta.bootstrap(x, test, 1.5, 0.1, 500)
    res2 <- get.theta.typed.bootstrap(x, 1,2, 1.5, 0.1, 500)


    #should have 95% CI of 0,1 and mean/median of 0.5
    expect_that(as.numeric(quantile(res[,1], probs=c(.025,.975), na.rm=T)),
                equals(c(0,Inf)))
    expect_that(as.numeric(quantile(res2[,1], probs=c(.025,.975), na.rm=T)),
                equals(c(0,Inf)))


})


test_that("performs correctly for test case 2 (points on a line)", {

    x<-rbind(c(1,0,0), c(2,1,0), c(2,-1,0), c(3,2,0),
             c(2,-2,0), c(3,3,0),c(3,-3,0))

    colnames(x) <-c("type","x","y")

    test <- function(a,b) {
        if (a[1] != 1) return(3)
        if (b[1] == 2) return(1)
        return(2)
    }

    #the medians for the null distribution should be 1,0.5,0
    res <- get.theta.bootstrap(x, test, c(1.5,2.5,3.5), c(0,1.5,2.5), 500)
    res2 <- get.theta.typed.bootstrap(x, 1, 2, c(1.5,2.5,3.5), c(0,1.5,2.5), 500)

    expect_that(median(res[,1], na.rm=T), equals(Inf))
    expect_that(median(res[,2], na.rm=T), equals(1))
    expect_that(median(res[,3], na.rm=T), equals(0))

    expect_that(median(res2[,1], na.rm=T), equals(Inf))
    expect_that(median(res2[,2], na.rm=T), equals(1))
    expect_that(median(res2[,3], na.rm=T), equals(0))


    #FIRST RANGE
    #deterministically Inf
    expect_that(as.numeric(quantile(res[,1], probs=c(.025,.975), na.rm=T)),
                equals(c(Inf,Inf)))
    expect_that(as.numeric(quantile(res2[,1], probs=c(.025,.975), na.rm=T)),
                equals(c(Inf,Inf)))

    #SECOND RANGE...should be 0 and Inf respectively a fairly large % of the time
    expect_that(as.numeric(quantile(res[,2], probs=c(0.025,.975), na.rm=T)),
                equals(c(0,Inf)))
    expect_that(as.numeric(quantile(res2[,2], probs=c(.025,.975), na.rm=T)),
                equals(c(0,Inf)))

    #THIRD RANGE
    #deterministically 0
    expect_that(as.numeric(quantile(res[,3], probs=c(.025,.975), na.rm=T)),
                equals(c(0,0)))
    expect_that(as.numeric(quantile(res2[,3], probs=c(.025,.975), na.rm=T)),
                equals(c(0,0)))


})


test_that ("fails nicely if x and y column names are not provided", {

    x<-cbind(rep(c(1,2),500), a=runif(1000,0,100), b=runif(1000,0,100))

    test <- function(a,b) {
        if (a[1] != 2) return(3)
        if (b[1] == 3) return(1)
        return(2)
    }

    expect_that(get.theta.bootstrap(x,test,seq(10,50,10), seq(0,40,10),100),
                throws_error("unique x and y columns must be defined"))

    expect_that(get.theta.ci(x,test,seq(10,50,10), seq(0,40,10),100),
                throws_error("unique x and y columns must be defined"))
})
