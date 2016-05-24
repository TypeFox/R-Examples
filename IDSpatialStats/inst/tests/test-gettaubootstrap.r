context("get.tau.bootstrap")

test_that("get.tau.bootstrap runs and returs 1 when it should", {

    x<-cbind(rep(c(1,2),50), x=runif(100,0,100), y=runif(100,0,100))

    colnames(x) <-c("type","x","y")

    test <- function(a,b) {return(1)}

    ########### REPRESENTATIVE
    #should return a matrix of all ones
    res <- get.tau.bootstrap(x, test, seq(10,100,10), seq(0,90,10), 10)
    expect_that(sum(res!=1),equals(0))
    expect_that(nrow(res),equals(10))


    ########### INDEPENDENT
    res <- get.tau.bootstrap(x, test, seq(10,100,10), seq(0,90,10), 10,
                             comparison.type="independent")
    expect_that(sum(res!=1),equals(0))
    expect_that(nrow(res),equals(10))

})


test_that("performs correctly for test case 1 (equilateral triangle)", {
    x <- rbind(c(1,0,0), c(1,1,0),c(2,.5,sqrt(.75)))
    colnames(x) <-c("type","x","y")

    test <- function(a,b) {
        if (a[1] != 1) return(3)
        if (b[1] == 2) return(1)
        return(2)
    }

    ########### REPRESENTATIVE
    res <- get.tau.bootstrap(x, test, 1.5, 0.1, 500)
    res2 <- get.tau.typed.bootstrap(x, 1,2, 1.5, 0.1, 500)

    #should have 95% CI of 1,1
    expect_that(as.numeric(quantile(res[,1], probs=c(.025,.975), na.rm=T)),
                equals(c(1,1)))

    expect_that(as.numeric(quantile(res2[,1], probs=c(.025,.975), na.rm=T)),
                equals(c(1,1)))

    ########### INDEPENDENT
    res <- get.tau.bootstrap(x, test, 1.5, 0.1, 500,
                             comparison.type="independent")
    res2 <- get.tau.typed.bootstrap(x, 1,2, 1.5, 0.1, 500,
                                    comparison.type="independent")

    #should have 95% CI of 1,1
    expect_that(as.numeric(quantile(res[,1], probs=c(.025,.975), na.rm=T)),
                equals(c(1,1)))

    expect_that(as.numeric(quantile(res2[,1], probs=c(.025,.975), na.rm=T)),
                equals(c(1,1)))

})


test_that("performs correctly for test case 2 (points on a line) - representative comparison group", {

    x<-rbind(c(1,0,0), c(2,1,0), c(2,-1,0), c(3,2,0),
             c(2,-2,0), c(3,3,0),c(3,-3,0))

    colnames(x) <-c("type","x","y")

    test <- function(a,b) {
        if (a[1] != 1) return(3)
        if (b[1] == 2) return(1)
        return(2)
    }

    ########### REPRESENTATIVE
    #the medians for the null distribution should be 2,1,0
    res <- get.tau.bootstrap(x, test, c(1.5,2.5,3.5), c(0,1.5,2.5), 1500)
    res2 <- get.tau.typed.bootstrap(x, 1, 2, c(1.5,2.5,3.5), c(0,1.5,2.5), 1500)

    expect_that(median(res[,1], na.rm=T), equals(2))
    expect_that(median(res[,2], na.rm=T), equals(1))
    expect_that(median(res[,3], na.rm=T), equals(0))

    expect_that(median(res2[,1], na.rm=T), equals(2))
    expect_that(median(res2[,2], na.rm=T), equals(1))
    expect_that(median(res2[,3], na.rm=T), equals(0))




    #FIRST RANGE
    #max would be only 1 type 2 used and in range = 1/(1/6) = 6...should occur
    #more than 2.5% of time
    #min would be 1, occuring just over .01% of the time
    expect_that(as.numeric(quantile(res[,1], probs=c(.001,.975), na.rm=T)),
                equals(c(1,6)))
    expect_that(as.numeric(quantile(res2[,1], probs=c(.001,.975), na.rm=T)),
                equals(c(1,6)))

    #SECOND RANGE
    #max would be 6, should occur less than 1% of the time
    #min should be 0, should occur 2.5% of the time
    expect_that(as.numeric(quantile(res[,2], probs=c(.025), na.rm=T)),
                equals(0))
    expect_that(as.numeric(quantile(res2[,2], probs=c(.025), na.rm=T)),
                equals(0))

    expect_that(as.numeric(quantile(res[,2], probs=c(.99), na.rm=T))<6,
                is_true())
    expect_that(as.numeric(quantile(res2[,2], probs=c(.99), na.rm=T))<6,
                is_true())



    #THIRD RANGE
    #Should be determinsitically 0 or NaN
    expect_that(as.numeric(quantile(res[,3], probs=c(.025,.975), na.rm=T)),
                equals(c(0,0)))
    expect_that(as.numeric(quantile(res2[,3], probs=c(.025,.975), na.rm=T)),
                equals(c(0,0)))


})

test_that("performs correctly for test case 2 (points on a line) - independent comparison group", {

    x<-rbind(c(1,0,0), c(2,1,0), c(2,-1,0), c(3,2,0),
             c(2,-2,0), c(3,3,0),c(3,-3,0))

    colnames(x) <-c("type","x","y")

    test <- function(a,b) {
        if (a[1] != 1) return(3)
        if (b[1] == 2) return(1)
        return(2)
    }

    ########### INDEPENDENT
    #the medians for the null distribution should be Inf,1,0
    res <- get.tau.bootstrap(x, test, c(1.5,2.5,3.5), c(0,1.5,2.5), 1500,
                             comparison.type="independent")
    res2 <- get.tau.typed.bootstrap(x, 1, 2, c(1.5,2.5,3.5), c(0,1.5,2.5), 1500,
                                    comparison.type="independent")

    expect_that(median(res[,1], na.rm=T), equals(Inf))
    expect_that(median(res[,2], na.rm=T), equals(1))
    expect_that(median(res[,3], na.rm=T), equals(0))

    expect_that(median(res2[,1], na.rm=T), equals(Inf))
    expect_that(median(res2[,2], na.rm=T), equals(1))
    expect_that(median(res2[,3], na.rm=T), equals(0))




    #FIRST RANGE
    #max would be Inf, occuring most of the time
    #min would be 1, occuring just over .01% of the time
    expect_that(as.numeric(quantile(res[,1], probs=c(.001,.975), na.rm=T)),
                equals(c(1,Inf)))
    expect_that(as.numeric(quantile(res2[,1], probs=c(.001,.975), na.rm=T)),
                equals(c(1,Inf)))

    #SECOND RANGE
    #max would be Inf, should occur around 25% of the time. .7 should be
    # reliably less than
    #min should be 0, should occur 2.5% of the time
    expect_that(as.numeric(quantile(res[,2], probs=c(.025), na.rm=T)),
                equals(0))
    expect_that(as.numeric(quantile(res2[,2], probs=c(.025), na.rm=T)),
                equals(0))

    expect_that(as.numeric(quantile(res[,2], probs=c(.7), na.rm=T))!=Inf,
                is_true())
    expect_that(as.numeric(quantile(res2[,2], probs=c(.7), na.rm=T))!=Inf,
                is_true())



    #THIRD RANGE
    #Should be determinsitically 0 or NaN
    expect_that(as.numeric(quantile(res[,3], probs=c(.025,.975), na.rm=T)),
                equals(c(0,0)))
    expect_that(as.numeric(quantile(res2[,3], probs=c(.025,.975), na.rm=T)),
                equals(c(0,0)))


})


test_that("get.tau.ci returns bootstrap cis when same seed", {
    x<-cbind(rep(c(1,2),50), x=runif(100,0,100), y=runif(100,0,100))

    colnames(x) <-c("type","x","y")

    test <- function(a,b) {
        if (a[1] != 1) return(3)
        if (b[1] == 2) return(1)
        return(2)
    }

    ####REPRESENTATIVE
    set.seed(787)
    res <- get.tau.bootstrap(x, test, seq(15,45,15), seq(0,30,15), 20)

    set.seed(787)
    ci1 <- get.tau.ci(x, test, seq(15,45,15), seq(0,30,15), 20)

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

    ### INDEPENDENT
    set.seed(787)
    res <- get.tau.bootstrap(x, test, seq(15,45,15), seq(0,30,15), 20,
                             comparison.type="independent")

    set.seed(787)
    ci1 <- get.tau.ci(x, test, seq(15,45,15), seq(0,30,15), 20,
                      comparison.type="independent")

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


test_that("fails nicely if x and y column names are not provided", {

    x<-cbind(rep(c(1,2),500), a=runif(1000,0,100), b=runif(1000,0,100))

    test <- function(a,b) {
        if (a[1] != 2) return(3)
        if (b[1] == 3) return(1)
        return(2)
    }

    expect_that(get.tau.bootstrap(x,test,seq(10,50,10), seq(0,40,10),100),
                throws_error("unique x and y columns must be defined"))

    expect_that(get.tau.ci(x,test,seq(10,50,10), seq(0,40,10),100),
                throws_error("unique x and y columns must be defined"))
})


##################DEPRECATED TESTS...TAKE TO LONG...NOW USING SMALLER CANONICAL
##################TESTS THAT HAVE VALUES THAT CAN BE WORKED OUT BY HAND


## test_that("CIs calculated from get.tau.bootstrap include the true value", {
##     set.seed(777)

##     x<-cbind(rep(c(1,2),250), x=runif(500,0,100), y=runif(500,0,100))

##     colnames(x) <-c("type","x","y")

##     test <- function(a,b) {
##         if (a[1] != 1) return(3)
##         if (b[1] == 1) return(1)
##         return(2)
##     }

##     res <- get.tau.ci(x, test, seq(10,100,10), seq(0,90,10), 200)

##     #print(res)

##     expect_that(sum(!(res[1,]<res[2,])),equals(0))
##     expect_that(sum(!(res[1,]<1)),equals(0))
##     expect_that(sum(!(res[2,]>1)),equals(0))

##     #repeat for typed data
##     res <- get.tau.typed.bootstrap(x, typeA=1, typeB=1,
##                                    seq(10,100,10), seq(0,90,10), 200)

##     ci <- matrix(nrow=2, ncol=ncol(res))

##     for (i in 1:ncol(ci)) {
##         ci[,i] <- quantile(res[,i], probs=c(0.025, 0.975))
##     }

##     res <- ci

##     expect_that(sum(!(res[1,]<res[2,])),equals(0))
##     expect_that(sum(!(res[1,]<1)),equals(0))
##     expect_that(sum(!(res[2,]>1)),equals(0))

## })

