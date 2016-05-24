context("get.theta.permute")

test_that("get.theta.permute returns appropriate values for test case 1 (equilateral triangle)" ,{

    x <- rbind(c(1,0,0), c(1,1,0),c(2,.5,sqrt(.75)))
    colnames(x) <-c("type","x","y")

    test <- function(a,b) {
        if (a[1] != 1) return(3)
        if (b[1] == 2) return(1)
        return(2)
    }


    #should return 1 for every permutation
    res <- get.theta.permute(x, test, 1.5, 0, 500)
    res2 <- get.theta.typed.permute(x, 1, 2, 1.5, 0, 500)


    expect_that(as.numeric(res), equals(rep(1,500)))
    expect_that(as.numeric(res2), equals(rep(1,500)))

})



test_that("get.theta.permute returns appropriate values for test case 2 (points on a line)" ,{
    x<-rbind(c(1,0,0), c(2,1,0), c(2,-1,0), c(3,2,0),
             c(2,-2,0), c(3,3,0),c(3,-3,0))

    colnames(x) <-c("type","x","y")

    test <- function(a,b) {
        if (a[1] != 1) return(3)
        if (b[1] == 2) return(1)
        return(2)
    }

    #the median of the null distribution should be 1 (includes infs so
    #  mean does not work)
    #the 95% CI equals 0,Inf with windows
    res <- get.theta.permute(x, test, c(1.5,2.5,3.5), c(0,1.5,2.5), 500)
    res2 <- get.theta.typed.permute(x, 1, 2, c(1.5,2.5,3.5), c(0,1.5,2.5), 500)

    expect_that(apply(res,2,median,na.rm=T), equals(rep(1,3), tolerance=0.1))
    expect_that(apply(res2, 2, median, na.rm=T), equals(rep(1,3), tolerance=0.1))

    for (i in 1:3) {
        expect_that(as.numeric(quantile(res[,i], probs=c(.025,.975))),
                    equals(c(0,Inf)))
        expect_that(as.numeric(quantile(res2[,i], probs=c(.025,.975))),
                    equals(c(0,Inf)))
    }

    #without windows the 95% CI should be 1/3 and 3
    res <- get.theta.permute(x, test, 4,0, 500)
    res2 <- get.theta.typed.permute(x, 1, 2, 4,0, 500)
    expect_that(as.numeric(quantile(res[,1], probs=c(.025,.975))),
                equals(c(1/3,3)))
    expect_that(as.numeric(quantile(res2[,1], probs=c(.025,.975))),
                equals(c(1/3,3)))
})


