context("Survival analysis related")

test_that("Survival data manipulation", {
    require(survival)

    failure_times <- rweibull(50, 2, 18)
    censoring_times <- 40*runif(50)
    events <- factor(ifelse(failure_times < censoring_times, 1, sample(c(0, 2:3), 50, TRUE)),
                     levels=0:3, labels=c("censored", "death", "competing event 1", "competing event 2"))
    ys <- Surv(pmin(failure_times, censoring_times), events)

    expect_equal(dichotomize(ys, to_factor=FALSE), ys[,"status"])
    expect_true(all(
        sign(table(dichotomize(ys), useNA="always") -
             table(dichotomize(ys, 30), useNA="always")) *
        c(1,1,1,1,-1) >= 0))

    # Missing value handling
    ys[1:2, "time"] <- NA
    ys[c(1, 3), "status"] <- NA
    expect_true(all(is.na(dichotomize(ys[1:3]))))
    expect_equal(dichotomize(ys[1:3], time=0, to_factor=FALSE), c(NA, NA, 0))
})
