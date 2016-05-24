test_that('aggregate', {
    partition <- matrix(c(0:10/10, 10:0/10, rep(0, 5),
                          rep(0, 5), 0:10/10, 10:0/10,
                          0:12/12, 1, 12:0/12),
                        byrow=FALSE,
                        ncol=3)
    colnames(partition) <- c('a', 'b', 'c')

    expect_equal(aggregate(c('a', 'b'),
                           c(1, 1),
                           partition),
                 pmin(partition[, 1], partition[, 2]))

    expect_equal(aggregate(c('a', 'b'),
                           c(0, 0),
                           partition),
                 rep(1, 27))

    expect_equal(aggregate(c('a', 'b'),
                           c(0.5, 0.8),
                           partition),
                 pmin(1, partition[, 1] + 0.5, partition[, 2] + 0.2))
})
