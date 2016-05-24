context("Unit tests of the pooledX functions")

for (i in 1:2) {
  if (i == 1) {
    pooledX <- pooledS

    # Simulate 4 covariance matrices
    n <- c(10, 4, 5, 7)
    sl <- createS(n, p = 7)
  } else {
    pooledX <- pooledP

    # Simulate 4 precision matrices
    n <- c(10, 8, 9, 13)
    sl <- lapply(createS(n, p = 7), solve)
  }

  test_that(sprintf("pooled%s works as intended", switch(i,"S","P")), {

    res <- pooledX(sl, n)

    expect_that(res, is_a("matrix"))
    expect_that(dim(res), equals(dim(sl[[1]])))
    expect_that(dimnames(res), equals(dimnames(sl[[1]])))

    # Length 1 argument
    expect_that(pooledX(sl[1], n[1]), equals(sl[[1]]))
  })

  test_that(sprintf("pooled%s's mle argument works as intended",
                    switch(i,"S","P")), {

    res1 <- pooledX(sl, n, mle = TRUE)
    res2 <- pooledX(sl, n, mle = FALSE)
    if (i == 1) {
      man1 <- Reduce(`+`, mapply(`*`, sl, n, SIMPLIFY = FALSE))/sum(n)
      man2 <- Reduce(`+`, mapply(`*`, sl, n-1, SIMPLIFY = FALSE))/sum(n-1)
    } else {
      tmp <- lapply(sl, solve)
      man1 <- solve(Reduce(`+`, mapply(`*`, tmp, n, SIMPLIFY = FALSE))/sum(n))
      man2 <-
        solve(Reduce(`+`, mapply(`*`, tmp, n-1, SIMPLIFY = FALSE))/sum(n-1))
    }

    # Standard
    expect_that(res1, is_a("matrix"))
    expect_that(dim(res1), equals(dim(sl[[1]])))
    expect_that(dimnames(res1), equals(dimnames(sl[[1]])))

    # Check equality
    expect_that(res1, equals(man1))
    expect_that(res2, equals(man2))

    # Check non-standard entries
    expect_error(pooledX(sl, n, mle = "A"))

  })

  test_that(sprintf("pooled%s's subset argument works as intended",
                    switch(i,"S","P")), {

    subset <- sample(c(TRUE, FALSE, FALSE, TRUE))
    res <- pooledX(sl, n, subset = subset)
    man <- pooledX(sl[subset], n[subset])

    # Standard
    expect_that(res, is_a("matrix"))
    expect_that(dim(res), equals(dim(sl[[1]])))
    expect_that(dimnames(res), equals(dimnames(sl[[1]])))

    # Check equality
    expect_that(res, equals(man))

  })

}
