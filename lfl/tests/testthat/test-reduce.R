test_that('reduce', {
    x <- matrix(c(1,1,1,1,0,0,1,1,0,0,
                  0,0,1,1,1,0,0,0,1,1,
                  1,1,1,1,1,1,0,0,0,0),
                byrow=FALSE, ncol=3)
    colnames(x) <- c('a', 'b', 'c')
    
    v <- colnames(x)
    names(v) <- colnames(x)

    s <- matrix(0, nrow=3, ncol=3)
    colnames(s) <- colnames(x)
    rownames(s) <- colnames(x)

    d <- fsets(x, v, s)



    rules <- farules(rules=list(c('c', 'a'),
                                c('c', 'b')),
                     statistics=matrix(as.numeric(1:2), nrow=2, ncol=1))

    expect_equal(reduce(d, rules, 0.66),
                 farules(rules=list(c('c', 'a')),
                         statistics=matrix(as.numeric(1), nrow=1, ncol=1)))

    expect_equal(reduce(d, rules, 0.67),
                 farules(rules=list(c('c', 'a'),
                                    c('c', 'b')),
                         statistics=matrix(as.numeric(1:2), nrow=2, ncol=1)))
})
