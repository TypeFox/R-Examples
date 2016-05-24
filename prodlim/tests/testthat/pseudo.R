context("Construction of pseudovalues")
test_that("pseudo",{
    library(prodlim)
    library(pseudo)
    # comparison to pseudoci
    # make sure we get the same
    # results with both packages
    set.seed(17)
    N <- 200
    ddd <- SimCompRisk(200)
    ttt <- c(3,5,10)
    # ttt <- ddd$time
    fff <- prodlim(Hist(time,event)~1,data=ddd)
    system.time(jack <- with(ddd,pseudoci(time,event,ttt)))
    system.time({jack2 <- jackknife.competing.risks(fff,times=ttt)})
    ## check individual 2
    expect_true(all(round(jack2[,2],9)==round(jack[[3]]$cause1[,2],9)))
    ## check all individuals
    expect_true(all(sapply(1:N,function(x){
        a <- round(jack[[3]]$cause1[x,],8)
        b <- round(jack2[x,],8)
        # all(a[!is.na(a)]==b[!is.na(b)])
        all(a[!is.na(a)]==b[!is.na(a)])
    })))
    ## the pseudoci function seems only slightly slower
    ## for small sample sizes (up to ca. 200) but
    ## much slower for large sample sizes:
    set.seed(17)
    N <- 200
    ddd <- SimCompRisk(200)
    ttt <- c(3,5,10)
    # ttt <- ddd$time
    fff <- prodlim(Hist(time,event)~1,data=ddd)
    system.time(jack <- with(ddd,pseudoci(time,event,ttt)))
    system.time({jack2 <- jackknife.competing.risks(fff,times=ttt)})
    expect_true(all(round(jack2[,1],9)==round(jack$pseudo$cause1[,1],9)))
    set.seed(17)
    N <- 2000
    ddd <- SimCompRisk(2000)
    ttt <- c(3,5,10)
    fff <- prodlim(Hist(time,event)~1,data=ddd)
    a <- system.time(jack <- with(ddd,pseudoci(time,event,ttt)))
    b <- system.time({jack2 <- jackknife.competing.risks(fff,times=ttt)})
    expect_less_than(a,b)
})

