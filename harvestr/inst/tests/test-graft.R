library(MCMCpack)
library(dostats)
library(plyr)

context("Grafting")
test_that("graft is reproducible", {
    results <- 
    replicate(2, simplify=F, {
        datasets <- farm(gather(3, seed = 20120604), {
            x1 <- rnorm(100)
            x2 <- rnorm(100)
            y <- rbinom(100, 1, p=plogis(x1 + x2))
            data.frame(y, x1, x2)
        })

        substreams <- llply(datasets, graft, 10)
        subchains <- harvest(substreams[[1]], MCMCregress
                            , formula=y~x1+x2, n=100)
    })
    expect_identical(noattr(results[[1]]), noattr(results[[2]]))
})

