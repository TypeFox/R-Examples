test_that('basic goedel composition', {
    R <- matrix(c(0.1, 0.6, 1, 0, 0, 0,
                  0, 0.3, 0.7, 0.9, 1, 1,
                  0, 0, 0.6, 0.8, 1, 0,
                  0, 1, 0.5, 0, 0, 0,
                  0, 0, 1, 1, 0, 0), byrow=TRUE, nrow=5)

    S <- matrix(c(0.9, 1, 0.9, 1,
                  1, 1, 1, 1,
                  0.1, 0.2, 0, 0.2,
                  0, 0, 0, 0,
                  0.7, 0.6, 0.5, 0.4,
                  1, 0.9, 0.7, 0.6), byrow=TRUE, nrow=6)

    RS <- matrix(c(0.6, 0.6, 0.6, 0.6,
                   1, 0.9, 0.7, 0.6,
                   0.7, 0.6, 0.5, 0.4,
                   1, 1, 1, 1,
                   0.1, 0.2, 0, 0.2), byrow=TRUE, nrow=5)
        
    expect_that(compose(R, S, alg='goedel', type='basic'), equals(RS))
})

test_that('basic goguen composition', {
    R <- matrix(c(0.1, 0.6, 1, 0, 0, 0,
                  0, 0.3, 0.7, 0.9, 1, 1,
                  0, 0, 0.6, 0.8, 1, 0,
                  0, 1, 0.5, 0, 0, 0,
                  0, 0, 1, 1, 0, 0), byrow=TRUE, nrow=5)

    S <- matrix(c(0.9, 1, 0.9, 1,
                  1, 1, 1, 1,
                  0.1, 0.2, 0, 0.2,
                  0, 0, 0, 0,
                  0.7, 0.6, 0.5, 0.4,
                  1, 0.9, 0.7, 0.6), byrow=TRUE, nrow=6)

    RS <- matrix(c(0.6, 0.6, 0.6, 0.6,
                   1, 0.9, 0.7, 0.6,
                   0.7, 0.6, 0.5, 0.4,
                   1, 1, 1, 1,
                   0.1, 0.2, 0, 0.2), byrow=TRUE, nrow=5)
        
    expect_that(compose(R, S, alg='goguen', type='basic'), equals(RS))
})

test_that('basic lukasiewicz composition', {
    R <- matrix(c(0.4, 0.9, 1, 0.9, 1, 0.6,
                  0.9, 0.5, 0, 0, 0, 0.2,
                  0, 0.9, 0.5, 1, 0.9, 1,
                  0.2, 1, 0, 0.2, 0.5, 0,
                  0.3, 1, 0.6, 1, 0.9, 0.9), byrow=TRUE, nrow=5)

    S <- matrix(c(0.1, 0.9, 0.2, 0,
                  0.8, 0.4, 1, 1,
                  0, 0, 0.8, 0.6,
                  0.4, 0, 1, 1,
                  0.9, 0.1, 1, 1,
                  0.2, 0, 0.6, 1), byrow=TRUE, nrow=6)

    RS <- matrix(c(0.9, 0.3, 1, 1,
                   0.3, 0.8, 0.5, 0.5,
                   0.8, 0.3, 1, 1,
                   0.8, 0.4, 1, 1,
                   0.8, 0.4, 1, 1), byrow=TRUE, nrow=5)

    expect_that(compose(R, S, alg='lukas', type='basic'), equals(RS))
})

test_that('lukasiewicz bandler-kohout subproduct composition', {
    R <- matrix(c(0.4, 0.9, 1, 0.9, 1, 0.6,
                  0.9, 0.5, 0, 0, 0, 0.2,
                  0, 0.9, 0.5, 1, 0.9, 1,
                  0.2, 1, 0, 0.2, 0.5, 0,
                  0.3, 1, 0.6, 1, 0.9, 0.9), byrow=TRUE, nrow=5)

    S <- matrix(c(0.1, 0.9, 0.2, 0,
                  0.8, 0.4, 1, 1,
                  0, 0, 0.8, 0.6,
                  0.4, 0, 1, 1,
                  0.9, 0.1, 1, 1,
                  0.2, 0, 0.6, 1), byrow=TRUE, nrow=6)

    RS <- matrix(c(0, 0, 0.8, 0.6,
                   0.2, 0.8, 0.3, 0.1,
                   0.2, 0, 0.6, 1,
                   0.8, 0.4, 1, 0.8,
                   0.3, 0, 0.7, 0.7), byrow=TRUE, nrow=5)

    expect_that(compose(R, S, alg='lukas', type='sub'), equals(RS))
})



test_that('lukasiewicz bandler-kohout superproduct composition', {
    R <- matrix(c(0.4, 0.9, 1, 0.9, 1, 0.6,
                  0.9, 0.5, 0, 0, 0, 0.2,
                  0, 0.9, 0.5, 1, 0.9, 1,
                  0.2, 1, 0, 0.2, 0.5, 0,
                  0.3, 1, 0.6, 1, 0.9, 0.9), byrow=TRUE, nrow=5)

    S <- matrix(c(0.1, 0.9, 0.2, 0,
                  0.8, 0.4, 1, 1,
                  0, 0, 0.8, 0.6,
                  0.4, 0, 1, 1,
                  0.9, 0.1, 1, 1,
                  0.2, 0, 0.6, 1), byrow=TRUE, nrow=6)

    RS <- matrix(c(1, 0.5, 0.9, 0.6,
                   0.1, 0.9, 0, 0,
                   0.9, 0.1, 0.7, 0.9,
                   0.6, 0.3, 0.2, 0,
                   1, 0.4, 0.8, 0.9), byrow=TRUE, nrow=5)

        
    expect_that(compose(R, S, alg='lukas', type='super'), equals(RS))
})



test_that('lukasiewicz bandler-kohout square composition', {
    R <- matrix(c(0.4, 0.9, 1, 0.9, 1, 0.6,
                  0.9, 0.5, 0, 0, 0, 0.2,
                  0, 0.9, 0.5, 1, 0.9, 1,
                  0.2, 1, 0, 0.2, 0.5, 0,
                  0.3, 1, 0.6, 1, 0.9, 0.9), byrow=TRUE, nrow=5)

    S <- matrix(c(0.1, 0.9, 0.2, 0,
                  0.8, 0.4, 1, 1,
                  0, 0, 0.8, 0.6,
                  0.4, 0, 1, 1,
                  0.9, 0.1, 1, 1,
                  0.2, 0, 0.6, 1), byrow=TRUE, nrow=6)

    RS <- matrix(c(0, 0, 0.8, 0.6,
                   0.1, 0.8, 0, 0,
                   0.2, 0, 0.6, 0.9,
                   0.6, 0.3, 0.2, 0,
                   0.3, 0, 0.7, 0.7), byrow=TRUE, nrow=5)

        
    expect_that(compose(R, S, alg='lukas', type='square'), equals(RS))
})

test_that('basic goedel composition with e', {
    R <- matrix(c(0.1, 0.6, 1, 0, 0, 0,
                  0, 0.3, 0.7, 0.9, 1, 1,
                  0, 0, 0.6, 0.8, 1, 0,
                  0, 1, 0.5, 0, 0, 0,
                  0, 0, 1, 1, 0, 0), byrow=TRUE, nrow=5)

    S <- matrix(c(0.9, 1, 0.9, 1,
                  1, 1, 1, 1,
                  0.1, 0.2, 0, 0.2,
                  0, 0, 0, 0,
                  0.7, 0.6, 0.5, 0.4,
                  1, 0.9, 0.7, 0.6), byrow=TRUE, nrow=6)

    E <- matrix(c(0, 0, 0, 0,
                  0, 0, 0, 0,
                  0, 0, 0, 0,
                  0, 0, 0, 0,
                  0, 0, 0, 0,
                  0, 0, 0, 0), byrow=TRUE, nrow=6)

    RS <- matrix(c(0.6, 0.6, 0.6, 0.6,
                   1, 0.9, 0.7, 0.6,
                   0.7, 0.6, 0.5, 0.4,
                   1, 1, 1, 1,
                   0.1, 0.2, 0, 0.2), byrow=TRUE, nrow=5)
        
    expect_that(compose(R, S, E, alg='goedel', type='basic'), equals(RS))
})


