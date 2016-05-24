context("Result extractors")

test_that("Parameter tuning", {
    # Simpler uses are covered in the examples section
    procedure <- list(
        RF1 = modeling_procedure("randomForest",
            parameter = list(mtry = c(1, 2, 3),
                             nodesize = c(1, 4, 10))),
        RF2 = modeling_procedure("randomForest",
            parameter = list(maxnodes = c(3, 5, 10)))
    )
    model <- fit(procedure, iris[-5], iris$Species, resample=cv)
    result <- evaluate(procedure, iris[-5], iris$Species, resample=cv,
                       .save=c(model=TRUE))
    tuning <- get_tuning(result)
    expect_is(tuning, "data.frame")
    expect_equal(sort(colnames(tuning)),
        c("error", "fold", "method", "parameter", "tuning_fold", "value"))
})

