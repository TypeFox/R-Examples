context("get.tau.permute")

test_that("get.tau.permute returns appropriate values for test case 1 (equilateral triangle)" ,{

    x <- rbind(c(1,0,0), c(1,1,0),c(2,.5,sqrt(.75)))
    colnames(x) <-c("type","x","y")

    test <- function(a,b) {
        if (a[1] != 1) return(3)
        if (b[1] == 2) return(1)
        return(2)
    }


    ###REPRESENTATIVE
    #should return 1 for every permutation
    res <- get.tau.permute(x, test, 1.5, 0, 500)
    res2 <- get.tau.typed.permute(x, 1, 2, 1.5, 0, 500)


    expect_that(as.numeric(res), equals(rep(1,500)))
    expect_that(as.numeric(res2), equals(rep(1,500)))

    ###INDEPENDENT
    #should return 1 for every permutation
    res <- get.tau.permute(x, test, 1.5, 0, 500,
                           comparison.type="independent")
    res2 <- get.tau.typed.permute(x, 1, 2, 1.5, 0, 500,
                                  comparison.type="independent")


    expect_that(as.numeric(res), equals(rep(1,500)))
    expect_that(as.numeric(res2), equals(rep(1,500)))

})


test_that("get.tau.permute returns appropriate values for test case 2 (points on a line)" ,{
    x<-rbind(c(1,0,0), c(2,1,0), c(2,-1,0), c(3,2,0),
             c(2,-2,0), c(3,3,0),c(3,-3,0))

    colnames(x) <-c("type","x","y")

    test <- function(a,b) {
        if (a[1] != 1) return(3)
        if (b[1] == 2) return(1)
        return(2)
    }

    ####REPRESENTATIVE
    #the mean of the null distribution should be 1
    #the 95% CI equals 0,2 with windows
    res <- get.tau.permute(x, test, c(1.5,2.5,3.5), c(0,1.5,2.5), 500)
    res2 <- get.tau.typed.permute(x, 1, 2, c(1.5,2.5,3.5), c(0,1.5,2.5), 500)

    expect_that(colMeans(res,na.rm=T), equals(rep(1,3), tolerance=0.1))
    expect_that(colMeans(res2, na.rm=T), equals(rep(1,3), tolerance=0.1))

    for (i in 1:3) {
        expect_that(as.numeric(quantile(res[,i], probs=c(.025,.975))),
                    equals(c(0,2)))
        expect_that(as.numeric(quantile(res2[,i], probs=c(.025,.975))),
                    equals(c(0,2)))
    }

    ####INDEPENDENT
    #the mean of the null distribution should be 1
    #the 95% CI equals 0,Inf with windows
    res <- get.tau.permute(x, test, c(1.5,2.5,3.5), c(0,1.5,2.5), 500,
                           comparison.type="independent")
    res2 <- get.tau.typed.permute(x, 1, 2, c(1.5,2.5,3.5), c(0,1.5,2.5), 500,
                                  comparison.type="independent")


    for (i in 1:3) {
        expect_that(as.numeric(quantile(res[,i], probs=c(.025,.5,.975))),
                    equals(c(0,1,Inf)))
        expect_that(as.numeric(quantile(res2[,i], probs=c(.025,.5,.975))),
                    equals(c(0,1,Inf)))
    }

})



test_that ("fails nicely if x and y column names are not provided", {
    x<-cbind(rep(c(1,2),500), a=runif(1000,0,100), b=runif(1000,0,100))

    test <- function(a,b) {
        if (a[1] != 2) return(3)
        if (b[1] == 3) return(1)
        return(2)
    }

    expect_that(get.tau.permute(x,test,seq(10,50,10), seq(0,40,10),100),
                throws_error("unique x and y columns must be defined"))

})


##################DEPRECATED TESTS...TAKE TO LONG...NOW USING SMALLER CANONICAL
##################TESTS THAT HAVE VALUES THAT CAN BE WORKED OUT BY HAND
## test_that("get.tau.permute cis enclose get.tau when no clustering exists",
##           {
##               set.seed(788)

##               x<-cbind(rep(c(1,2),250), x=runif(500,0,100), y=runif(500,0,100))

##               colnames(x) <-c("type","x","y")

##               test <- function(a,b) {
##                   if (a[1] != 1) return(3)
##                   if (b[1] == 1) return(1)
##                   return(2)
##               }

##               #plot(x[,"x"],x[,"y"], col=x[,"type"])

##               res <- get.tau.permute(x, test, seq(10,100,10), seq(0,90,10), 200)
##               res2 <- get.tau(x, test, seq(10,100,10), seq(0,90,10))

##               for (i in 1:10) {
##                   tmp <- quantile(res[,i], probs=c(0.025, .975), na.rm=T)
##                   #print(res2[i])
##                   #print(tmp)
##                   expect_that(res2[i]>=tmp[1], is_true())
##                   expect_that(res2[i]<=tmp[2], is_true())
##               }
##           })

## test_that("get.tau.permute cis do not enclose get.tau at extremes when no clustering exists",
##           {
##               set.seed(788)

##               #first generate 200 random uniform points
##               x<-cbind(1, x=runif(300,0,100), y=runif(300,0,100))
##               colnames(x) <-c("type","x","y")

##               #add a seed point
##               x<-rbind(x,c(2,50,50))

##               #generate 200 normally distibuted points around this
##               x<-rbind(x,cbind(3,rnorm(300,50,20),rnorm(300,50,20)))

##               test <- function(a,b) {
##                   if (a[1] != 2) return(3)
##                   if (b[1] == 3) return(1)
##                   return(2)
##               }

##               ## res <- get.tau.permute(x,test,seq(10,50,10), seq(0,40,10), 200)
##               res2 <- get.tau(x,test,seq(10,50,10), seq(0,40,10))


##               for (i in c(1,5)) {
##                   tmp <- quantile(res[,i], probs=c(0.025, .975), na.rm=T)
##                   #print(tmp)
##                   #print(res2[i])
##                   expect_that((res2[i]>=tmp[1]) & (res2[i]<=tmp[2]) ,
##                               is_false())
##               }



##               res <- get.tau.typed.permute(x,2,3,seq(10,50,10),
##                                            seq(0,40,10), 100)


##               for (i in c(1,5)) {
##                   tmp <- quantile(res[,i], probs=c(0.025, .975), na.rm=T)
##                   #print(res2[i])
##                   #print(tmp)
##                   expect_that((res2[i]>=tmp[1]) & (res2[i]<=tmp[2]) ,
##                               is_false())
##               }


##           })


