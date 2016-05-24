library(rankdist)
context("LogC")

faiC = function(fai,t.lst=NULL){
    if (is.null(t.lst))
        t.lst = t.gen(length(fai))
    d = length(fai) # d = t - 1
    K = matrix(rep(0,d^2),ncol = d, nrow = d)
    for ( i in 1:d){
        K = -1 * fai[i] * t.lst[[i]] + K
    }
    K = exp(K)
    K[upper.tri(K)] = 0
    ones = rep(1,d)
    denom = exp(sum(log(rowSums(K) + ones)))
    denom
}

t.gen = function(d){
    t.lst = list()
    t.lst[[d]] = matrix(rep(1:d,d),ncol = d, nrow = d,byrow=T)
    left.mask = matrix(rep(0,d^2),ncol = d, nrow = d)
    left.mask[2:d,1:(d-1)] = diag(rep(1,d-1))
    t.lst[[d]][upper.tri(left.mask)] = 0
    for ( i in 1:(d-1)){
        t.lst[[d-i]] = left.mask%*%t.lst[[d-i+1]]
        diag(t.lst[[d-i]]) = c(rep(0,i),1:(d-i))
        t.lst[[d-i]][upper.tri(left.mask)] = 0
    }
    t.lst
}
a = runif(n=10)
b = runif(n=10)
c = runif(n=10)
d = runif(n=10)

test_that("constant is correct",{
    expect_equal(LogC(a),log(faiC(a,t.gen(length(a)))))
    expect_equal(LogC(b),log(faiC(b,t.gen(length(b)))))
    expect_equal(LogC(c),log(faiC(c,t.gen(length(c)))))
    expect_equal(LogC(d),log(faiC(d,t.gen(length(d)))))
})

