set.seed(34523)

test_that('minimum (Goedel) t-norm', {
    tnorm <- .tnorms[['goedel']]
    expect_that(tnorm(0.2, 0.5, 0.1, 0.3), equals(0.1))
    expect_that(tnorm(0.4, 0.5, 0.3), equals(0.3))
    expect_that(tnorm(0.2, 0.5, 0.9), equals(0.2))
    expect_that(tnorm(0.2, 0.5, 0.0), equals(0))
    expect_that(tnorm(1, 1, 1, 1), equals(1))
    expect_that(tnorm(1, 0.9, 1, 1), equals(0.9))
    expect_that(tnorm(0.2, NA, 1, na.rm=TRUE), equals(0.2))
    expect_that(tnorm(0.2, NA, 0, na.rm=TRUE), equals(0))

    expect_that(tnorm(c(0.2, 0.5, 0.1, 0.3)), equals(0.1))
    expect_that(tnorm(c(0.4, 0.5, 0.3)), equals(0.3))
    expect_that(tnorm(c(0.2, 0.5, 0.9)), equals(0.2))
    expect_that(tnorm(c(0.2, 0.5, 0.0)), equals(0))
    expect_that(tnorm(c(1, 1, 1, 1)), equals(1))
    expect_that(tnorm(c(1, 0.9, 1, 1)), equals(0.9))
    expect_that(tnorm(c(0.2, NA, 1), na.rm=TRUE), equals(0.2))
    expect_that(tnorm(c(0.2, NA, 0), na.rm=TRUE), equals(0))
})


test_that('lukasiewicz t-norm', {
    tnorm <- .tnorms[['lukasiewicz']]
    expect_that(tnorm(0.2, 0.5, 0.1, 0.3), equals(0))
    expect_that(tnorm(0.8, 0.5, 0.9), equals(0.2))
    expect_that(tnorm(1, 1, 1, 1), equals(1))
    expect_that(tnorm(1, 0.9, 1, 1), equals(0.9))
    expect_that(tnorm(0.2, 0.5, 0.0), equals(0))
    expect_that(tnorm(0.2, 0.4, NA), equals(0))
    expect_that(tnorm(1, 0.9, 1, 1, NA, na.rm=TRUE), equals(0.9))

    expect_that(tnorm(c(0.2, 0.5, 0.1, 0.3)), equals(0))
    expect_that(tnorm(c(0.8, 0.5, 0.9)), equals(0.2))
    expect_that(tnorm(c(1, 1, 1, 1)), equals(1))
    expect_that(tnorm(c(1, 0.9, 1, 1)), equals(0.9))
    expect_that(tnorm(c(0.2, 0.5, 0.0)), equals(0))
    expect_that(tnorm(c(0.2, 0.4, NA)), equals(0))
    expect_that(tnorm(c(1, 0.9, 1, 1, NA), na.rm=TRUE), equals(0.9))
})


test_that('product (goguen) t-norm', {
    tnorm <- .tnorms[['goguen']]
    expect_that(tnorm(0.2, 0.5, 0.1, 0.3), equals(0.2 * 0.5 * 0.1 * 0.3))
    expect_that(tnorm(0.8, 0.5, 0.9), equals(0.8 * 0.5 * 0.9))
    expect_that(tnorm(1, 1, 1, 1), equals(1))
    expect_that(tnorm(1, 0.9, 1, 1), equals(0.9))
    expect_that(tnorm(0.2, 0.5, 0.0), equals(0))
    expect_that(tnorm(1, 0.9, 1, NA, 1, na.rm=TRUE), equals(0.9))

    expect_that(tnorm(c(0.2, 0.5, 0.1, 0.3)), equals(0.2 * 0.5 * 0.1 * 0.3))
    expect_that(tnorm(c(0.8, 0.5, 0.9)), equals(0.8 * 0.5 * 0.9))
    expect_that(tnorm(c(1, 1, 1, 1)), equals(1))
    expect_that(tnorm(c(1, 0.9, 1, 1)), equals(0.9))
    expect_that(tnorm(c(0.2, 0.5, 0.0)), equals(0))
    expect_that(tnorm(c(1, 0.9, 1, NA, 1), na.rm=TRUE), equals(0.9))
})


for (ttt in names(.tnorms)) {
    test_that(paste(ttt, 't-norm borders'), {
        tnorm <- .tnorms[[ttt]]

        expect_that(tnorm(), equals(1))
        expect_that(tnorm(0.2, NA, 1), equals(as.numeric(NA)))
        expect_that(tnorm(0.2, NA, 0), equals(0))
        expect_that(tnorm(0.2, Inf, 0), throws_error('argument out of range 0..1'))
        expect_that(tnorm(0.2, -Inf, 0), throws_error('argument out of range 0..1'))
        expect_that(tnorm(0.2, 3, 0), throws_error('argument out of range 0..1'))
        expect_that(tnorm(0.2, -3, 0), throws_error('argument out of range 0..1'))
        expect_that(tnorm(0.2, NaN, 0), throws_error('NaN argument'))

        expect_that(tnorm(c()), equals(1))
        expect_that(tnorm(c(0.2, NA, 1)), equals(as.numeric(NA)))
        expect_that(tnorm(c(0.2, NA, 0)), equals(0))
        expect_that(tnorm(c(0.2, Inf, 0)), throws_error('argument out of range 0..1'))
        expect_that(tnorm(c(0.2, -Inf, 0)), throws_error('argument out of range 0..1'))
        expect_that(tnorm(c(0.2, 3, 0)), throws_error('argument out of range 0..1'))
        expect_that(tnorm(c(0.2, -3, 0)), throws_error('argument out of range 0..1'))
        expect_that(tnorm(c(0.2, NaN, 0)), throws_error('NaN argument'))
    })
}


test_that('parallel minimum t-norm', {
    tnorm <- .ptnorms[['goedel']]
    expect_that(tnorm(c(0.2, 0.5, 0.4, 0.9, 1),
                    c(0.5, 0.9, 0.8, 1.0, 1),
                    c(0.6, 0.4, 0.0, 1.0, 1),
                    c(0.3, 0.7, 0.5, 1.0, 1)),
                equals(c(0.2, 0.4, 0.0, 0.9, 1)))

    expect_that(tnorm(0.2, 0.5), equals(0.2))
    expect_that(tnorm(0.2, 0.5, 0.0), equals(0))
    expect_that(tnorm(c(0.2, 0.5, 0.0)), equals(c(0.2, 0.5, 0.0)))
})


test_that('parallel lukasiewicz t-norm', {
    tnorm <- .ptnorms[['lukasiewicz']]
    expect_that(tnorm(c(0.2, 0.8, 0.4, 0.9, 1),
                    c(0.5, 0.9, 0.8, 1.0, 1),
                    c(0.6, 0.5, 0.0, 1.0, 1),
                    c(0.3, 0.9, 0.5, 1.0, 1)),
                equals(c(0, 0.1, 0.0, 0.9, 1)))

    expect_that(tnorm(0.7, 0.8, 0.6), equals(0.1))
    expect_that(tnorm(0.2, 0.5, 0.0), equals(0))
    expect_that(tnorm(c(0.2, 0.5, 0.0)), equals(c(0.2, 0.5, 0.0)))
})


test_that('parallel product t-norm', {
    tnorm <- .ptnorms[['goguen']]
    expect_that(tnorm(c(0.2, 0.5, 0.4, 0.9, 1),
                    c(0.5, 0.9, 0.8, 1.0, 1),
                    c(0.6, 0.4, 0.0, 1.0, 1),
                    c(0.3, 0.7, 0.5, 1.0, 1)),
                equals(c(0.2 * 0.5 * 0.6 * 0.3,
                        0.5 * 0.9 * 0.4 * 0.7,
                        0.0, 0.9, 1)))

    expect_that(tnorm(0.2, 0.5), equals(0.2 * 0.5))
    expect_that(tnorm(0.2, 0.5, 0.0), equals(0))
    expect_that(tnorm(c(0.2, 0.5, 0.0)), equals(c(0.2, 0.5, 0.0)))
})


for (ttt in names(.ptnorms)) {
    test_that(paste('parallel', ttt, 't-norm borders'), {
        tnorm <- .ptnorms[[ttt]]

        expect_true(is.null(tnorm()))
        expect_that(tnorm(c(0.2, NA, 1), c(0.8, 0.6, NA))[2:3], equals(as.numeric(c(NA, NA))))
        expect_that(tnorm(c(0.2, NA, 0), c(0.8, 0, NA))[2:3], equals(c(0, 0)))
        expect_that(tnorm(c(0.2, 0.9, 0), c(0.8, 0, Inf)), throws_error('argument out of range 0..1'))
        expect_that(tnorm(c(0.2, 0.9, 0), c(0.8, 0, -Inf)), throws_error('argument out of range 0..1'))
        expect_that(tnorm(c(0.2, 0.9, 0), c(0.8, 0, 3)), throws_error('argument out of range 0..1'))
        expect_that(tnorm(c(0.2, 0.9, 0), c(0.8, 0, -3)), throws_error('argument out of range 0..1'))
        expect_that(tnorm(c(0.2, 0.9, 0), c(0.8, 0, NaN)), throws_error('NaN argument'))

        mr <- matrix(runif(12), nrow=3, ncol=4)
        colnames(mr) <- LETTERS[1:4]
        rownames(mr) <- letters[1:3]

        m0 <- matrix(0, nrow=3, ncol=4)
        colnames(m0) <- colnames(mr)
        rownames(m0) <- rownames(mr)
        expect_that(tnorm(mr, 0), equals(m0))
        expect_that(tnorm(mr, m0), equals(m0))

        m1 <- matrix(1, nrow=3, ncol=4)
        expect_that(tnorm(mr, 1), equals(mr))
        expect_that(tnorm(mr, m1), equals(mr))

        mx <- matrix(tnorm(c(mr), c(mr)), nrow=3, ncol=4)
        colnames(mx) <- colnames(mr)
        rownames(mx) <- rownames(mr)
        expect_that(tnorm(mr, mr), equals(mx))
    })
}


test_that('goedel t-conorm', {
    tconorm <- .tconorms[['goedel']]
    expect_that(tconorm(0.2, 0.5, 0.1, 0.3), equals(0.5))
    expect_that(tconorm(0.4, 0.5, 0.8), equals(0.8))
    expect_that(tconorm(0.9, 0.5, 0.2), equals(0.9))
    expect_that(tconorm(0.2, 1, 0.0), equals(1))
    expect_that(tconorm(0, 0, 0, 0), equals(0))
    expect_that(tconorm(0.2, NA, 0.5, na.rm=TRUE), equals(0.5))
    expect_that(tconorm(0.2, NA, 1, na.rm=TRUE), equals(1))

    expect_that(tconorm(c(0.2, 0.5, 0.1, 0.3)), equals(0.5))
    expect_that(tconorm(c(0.4, 0.5, 0.8)), equals(0.8))
    expect_that(tconorm(c(0.9, 0.5, 0.2)), equals(0.9))
    expect_that(tconorm(c(0.2, 1, 0.0)), equals(1))
    expect_that(tconorm(c(0, 0, 0, 0)), equals(0))
    expect_that(tconorm(c(0.2, NA, 0.5), na.rm=TRUE), equals(0.5))
    expect_that(tconorm(c(0.2, NA, 1), na.rm=TRUE), equals(1))
})


test_that('lukasiewicz t-conorm', {
    tconorm <- .tconorms[['lukasiewicz']]
    expect_that(tconorm(0.2, 0.5, 0.1, 0.0), equals(0.8))
    expect_that(tconorm(0.4, 0.5, 0.8), equals(1))
    expect_that(tconorm(1, 1, 1), equals(1))
    expect_that(tconorm(0, 0, 0, 0), equals(0))
    expect_that(tconorm(0.2, NA, 0.5, na.rm=TRUE), equals(0.7))
    expect_that(tconorm(0.2, NA, 1, na.rm=TRUE), equals(1))

    expect_that(tconorm(c(0.2, 0.5, 0.1, 0.0)), equals(0.8))
    expect_that(tconorm(c(0.4, 0.5, 0.8)), equals(1))
    expect_that(tconorm(c(1, 1, 1)), equals(1))
    expect_that(tconorm(c(0, 0, 0, 0)), equals(0))
    expect_that(tconorm(c(0.2, NA, 0.5), na.rm=TRUE), equals(0.7))
    expect_that(tconorm(c(0.2, NA, 1), na.rm=TRUE), equals(1))
})


test_that('goguen t-conorm', {
    tconorm <- .tconorms[['goguen']]
    expect_that(tconorm(0.2, 0.5, 0.1, 0.3), equals(0.748))
    expect_that(tconorm(0.2, 1, 0.0), equals(1))
    expect_that(tconorm(0, 0, 0, 0), equals(0))
    expect_that(tconorm(0.2, NA, 0.5, na.rm=TRUE), equals(0.6))
    expect_that(tconorm(0.2, NA, 1, na.rm=TRUE), equals(1))

    expect_that(tconorm(c(0.2, 0.5, 0.1, 0.3)), equals(0.748))
    expect_that(tconorm(c(0.2, 1, 0.0)), equals(1))
    expect_that(tconorm(c(0, 0, 0, 0)), equals(0))
    expect_that(tconorm(c(0.2, NA, 0.5), na.rm=TRUE), equals(0.6))
    expect_that(tconorm(c(0.2, NA, 1), na.rm=TRUE), equals(1))
})


for (ttt in names(.tconorms)) {
    test_that(paste(ttt, 't-conorm borders'), {
        tconorm <- .tconorms[[ttt]]

        expect_that(tconorm(), equals(0))
        expect_that(tconorm(0.2, NA, 0), equals(as.numeric(NA)))
        expect_that(tconorm(0.2, NA, 1), equals(1))
        expect_that(tconorm(0.2, Inf, 0), throws_error('argument out of range 0..1'))
        expect_that(tconorm(0.2, -Inf, 0), throws_error('argument out of range 0..1'))
        expect_that(tconorm(0.2, 3, 0), throws_error('argument out of range 0..1'))
        expect_that(tconorm(0.2, -3, 0), throws_error('argument out of range 0..1'))
        expect_that(tconorm(0.2, NaN, 0), throws_error('NaN argument'))

        expect_that(tconorm(c()), equals(0))
        expect_that(tconorm(c(0.2, NA, 0)), equals(as.numeric(NA)))
        expect_that(tconorm(c(0.2, NA, 1)), equals(1))
        expect_that(tconorm(c(0.2, Inf, 0)), throws_error('argument out of range 0..1'))
        expect_that(tconorm(c(0.2, -Inf, 0)), throws_error('argument out of range 0..1'))
        expect_that(tconorm(c(0.2, 3, 0)), throws_error('argument out of range 0..1'))
        expect_that(tconorm(c(0.2, -3, 0)), throws_error('argument out of range 0..1'))
        expect_that(tconorm(c(0.2, NaN, 0)), throws_error('NaN argument'))
    })
}


test_that('goedel residuum', {
    resid <- .residua[['goedel']]
    expect_that(resid(c(0, 0.2, 0.8, 1), 1), equals(c(1, 1, 1, 1)))
    expect_that(resid(c(0, 0.2, 0.8, 1), 0), equals(c(1, 0, 0, 0)))
    expect_that(resid(c(0, 0.2, 0.8, 1), 0.5), equals(c(1, 1, 0.5, 0.5)))
    expect_that(resid(c(0, 0.2, 0.8, 1), c(0.3, 0.9)), equals(c(1, 1, 0.3, 0.9)))
})



test_that('lukasiewicz residuum', {
    resid <- .residua[['lukasiewicz']]
    expect_that(resid(c(0, 0.2, 0.8, 1), 1), equals(c(1, 1, 1, 1)))
    expect_that(resid(c(0, 0.2, 0.8, 1), 0), equals(c(1, 0.8, 0.2, 0)))
    expect_that(resid(c(0, 0.2, 0.8, 1), 0.5), equals(c(1, 1, 0.7, 0.5)))
    expect_that(resid(c(0, 0.2, 0.8, 1), c(0.3, 0.9)), equals(c(1, 1, 0.5, 0.9)))
})


test_that('goguen residuum', {
    resid <- .residua[['goguen']]
    expect_that(resid(c(0, 0.2, 0.8, 1), 1), equals(c(1, 1, 1, 1)))
    expect_that(resid(c(0, 0.2, 0.8, 1), 0), equals(c(1, 0, 0, 0)))
    expect_that(resid(c(0, 0.2, 0.8, 1), 0.5), equals(c(1, 1, 0.625, 0.5)))
    expect_that(resid(c(0, 0.2, 0.8, 1), c(0.3, 0.9)), equals(c(1, 1, 0.375, 0.9)))
})


for (ttt in names(.residua)) {
    test_that(paste(ttt, 'residua borders'), {
        resid <- .residua[[ttt]]

        expect_that(resid(0.2, NA), equals(as.numeric(NA)))
        expect_that(resid(0, NA), equals(1))
        expect_that(resid(0.2, Inf), throws_error('argument out of range 0..1'))
        expect_that(resid(0.2, -Inf), throws_error('argument out of range 0..1'))
        expect_that(resid(0.2, 3), throws_error('argument out of range 0..1'))
        expect_that(resid(0.2, -3), throws_error('argument out of range 0..1'))
        expect_that(resid(0.2, NaN), throws_error('NaN argument'))

        expect_that(resid(NA, 0.2), equals(as.numeric(NA)))
        expect_that(resid(Inf, 0.2), throws_error('argument out of range 0..1'))
        expect_that(resid(-Inf, 0.2), throws_error('argument out of range 0..1'))
        expect_that(resid(3, 0.2), throws_error('argument out of range 0..1'))
        expect_that(resid(-3, 0.2), throws_error('argument out of range 0..1'))
        expect_that(resid(NaN, 0.2), throws_error('NaN argument'))
    })
}


test_that('involutive negation', {
    f <- .negations[['involutive']]
    expect_that(f(c(0, 0.2, NA, 0.8, 1)), equals(c(1, 0.8, NA, 0.2, 0)))

    m <- matrix(c(0, 0.2, NA, 0.8, 1, 0.3), nrow=2)
    colnames(m) <- letters[1:3]
    rownames(m) <- letters[1:2]
    r <- matrix(c(1, 0.8, NA, 0.2, 0, 0.7), nrow=2)
    colnames(r) <- letters[1:3]
    rownames(r) <- letters[1:2]
    expect_that(f(m), equals(r))

    expect_that(f(c(-3, 0.2)), throws_error('argument out of range 0..1'))
    expect_that(f(c(3, 0.2)), throws_error('argument out of range 0..1'))
    expect_that(f(c(-Inf, 0.2)), throws_error('argument out of range 0..1'))
    expect_that(f(c(Inf, 0.2)), throws_error('argument out of range 0..1'))
    expect_that(f(c(NaN, 0.2)), throws_error('NaN argument'))
})


test_that('strict negation', {
    f <- .negations[['strict']]
    expect_that(f(c(0, 0.2, NA, 0.8, 1)), equals(c(1, 0, NA, 0, 0)))

    m <- matrix(c(0, 0.2, NA, 0.8, 1, 0.3), nrow=2)
    colnames(m) <- letters[1:3]
    rownames(m) <- letters[1:2]
    r <- matrix(c(1, 0, NA, 0, 0, 0), nrow=2)
    colnames(r) <- letters[1:3]
    rownames(r) <- letters[1:2]
    expect_that(f(m), equals(r))

    expect_that(f(c(-3, 0.2)), throws_error('argument out of range 0..1'))
    expect_that(f(c(3, 0.2)), throws_error('argument out of range 0..1'))
    expect_that(f(c(-Inf, 0.2)), throws_error('argument out of range 0..1'))
    expect_that(f(c(Inf, 0.2)), throws_error('argument out of range 0..1'))
    expect_that(f(c(NaN, 0.2)), throws_error('NaN argument'))
})


