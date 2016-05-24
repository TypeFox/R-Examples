test_that('sel.farules', {
    orig <- farules(rules=list(letters[1:3],
                            letters[2:5],
                            letters[4],
                            letters[3:8]),
                    statistics=matrix(runif(16), nrow=4))

    res <- sel(orig, 2:3)
    expect_equal(res$rules, orig$rules[2:3])
    expect_equal(res$statistics, orig$statistics[2:3, ])

    res <- sel(orig, 1)
    expect_equal(res$rules, orig$rules[1])
    expect_equal(res$statistics, orig$statistics[1, , drop=FALSE])
})
