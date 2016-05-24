library(docopulae)
context('base functions')


## lproduct
test_that('lproduct', {
    expect_equal(lproduct(list()),
                 list())
    expect_equal(lproduct(list(list(1))),
                 list(list(1)))
    expect_equal(lproduct(list(list(1), list(2, 3))),
                 list(list(1, 2), list(1, 3)))
    expect_equal(lproduct(list(list(1, 2), list(3))),
                 list(list(1, 3), list(2, 3)))
    expect_equal(lproduct(list(list(list(1)))), # leave inner list(1) untouched
                 list(list(list(1))))
})


## clusterPeak
test_that('clusterPeak handles distance correctly', {
    expect_equal(clusterPeak(matrix(c(1, 2)), 2:1, 1),
                 c(1, 1))
    expect_equal(clusterPeak(matrix(c(1, 2 + 1e-15)), 2:1, 1),
                 c(1, 2))
})


test_that('clusterPeak works for a simple random example', {
    n = 100
    x1 = matrix(runif(n*2), ncol=2)
    x2 = matrix(runif(n*2), ncol=2)
    x3 = matrix(runif(n*2), ncol=2)
    x4 = matrix(runif(n*2), ncol=2)
    off = 1 + sqrt(2) + 1e-15
    x2 = sweep(x2, 2, c(off, 0), '+')
    x3 = sweep(x3, 2, c(0, off), '+')
    x4 = sweep(x4, 2, c(off, off), '+')
    x = rbind(x1, x2, x3, x4)
    tcl = rep(1:4, each=n) # true class
    ord = sample(1:nrow(x), nrow(x))
    x = x[ord,]
    tcl = tcl[ord]
    #plot(x[,1], x[,2], col=tcl)

    expect_true({
        cl = clusterPeak(x, nrow(x):1, sqrt(2))
        clg = split(cl, tcl) # groups
        all(sapply(clg, function(g) all(g == g[1]))) # are all groups consistent
    })
})


## roworder
test_that('roworder works in trivial cases', {
    expect_equal(roworder(matrix(nrow=0, ncol=1)), integer(0))
    expect_equal(roworder(data.frame(numeric(0))), integer(0))
    expect_equal(roworder(matrix(c(1, 2), ncol=1)), c(1, 2))
    expect_equal(roworder(data.frame(c(1, 2), ncol=1)), c(1, 2))
})


test_that('roworder works for expand.grid example', {
    x1 = 1:3
    x2 = 1:5
    x3 = 1:7
    x = as.matrix(expand.grid(x1, x2, x3))
    dimnames(x) = list()
    expect_equal(x[roworder(x),],
                 cbind(rep(x1, each=length(x2)*length(x3)),
                       rep(rep(x2, each=length(x3)), length(x1)),
                       rep(x3, length(x1)*length(x2))))
})


## rowmatch
test_that('rowmatch works in trivial cases', {
    a = matrix(nrow=0, ncol=2)
    b = matrix(1:4, ncol=2)
    expect_equal(rowmatch(a, b), integer(0))
    expect_equal(rowmatch(b, a), as.integer(c(NA, NA)))
    expect_equal(rowmatch(matrix(1), matrix(1)), 1)
    expect_equal(rowmatch(matrix(1:2), matrix(2:1)), 2:1)

    a = matrix(1:4, ncol=2, byrow=T)
    b = matrix(-1:6, ncol=2, byrow=T)
    expect_equal(rowmatch(a, b), 2:3)
    expect_equal(rowmatch(b, a), as.integer(c(NA, 1, 2, NA)))
})


test_that('rowmatch works for a random example', {
    ## subset
    x1 = matrix(rnorm(100*4), ncol=4)
    i = sample(1:nrow(x1), 50)
    expect_equal(rowmatch(x1[i,], x1), i)

    ## no match
    x2 = x1; x2[i,] = matrix(rnorm(50*4), ncol=4)
    expect_equal(which(is.na(rowmatch(x2, x1))), sort(i))
    expect_equal(c(na.omit(rowmatch(x2, x1))), (1:100)[-i])
})
    

### roworder and rowmatch
#test_that('rowmatch works using roworder and vice versa', {
    #x = matrix(rnorm(100*4), ncol=4)
    #i = sample(1:nrow(x), nrow(x) / 3)
    #y = x[i,]


#})


