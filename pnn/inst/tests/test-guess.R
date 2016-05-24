context("Guess")

test_that("Norms dataset", {
    a <- cbind(1,rnorm(10,1,1))
    b <- cbind(5,rnorm(10,5,1))
    trainset <- as.data.frame(rbind(a,b))
    pnn <- learn(set=trainset)
    pnn$sigma <- 0.3
    (res <- guess.probabilities.of.each.category(nn=pnn, X=1))
    expect_that(length(res), equals(2))
    expect_that(as.vector(res[1] > 0.5), is_true())
    expect_that(guess.category(nn=pnn, X=1), equals("1"))
})
