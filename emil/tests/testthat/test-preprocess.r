context("Pre-processing")

x <- sweep(matrix(0, 6, 6), 1, 1:6, "+")
na.ind <- arrayInd(sample(length(x)/2, 10), c(6,6))[,2:1]
x[na.ind] <- NA
y <- NULL
fold <- rep(0:1, each=3)

test_that("Median imputation", {
    sets <- pre_impute_median(pre_split(x, y, fold))
    expect_true(all(sets$test$x[na.ind] == 5))
})

test_that("k-NN imputation", {
    sets <- pre_impute_knn(pre_split(x, y, fold), k=1, distance_matrix=dist(x))
    expect_true(all(sets$test$x[na.ind] == 4))
    sets2 <- pre_impute_knn(pre_split(x, y, fold), k=1, distance_matrix="auto")
    expect_identical(sets, sets2)

    x[5,1:3] <- NA
    sets <- pre_impute_knn(pre_split(x, y, fold), k=2, distance_matrix="auto")
    expect_equivalent(sets$fit$x[2,], rep(5,6))
    expect_equivalent(sets$test$x[na.ind], ifelse(na.ind[,2] <= 3, 5, 4.5))
})

test_that("Factor to logical conversion", {
    df <- data.frame(
        email = factor(sample(2, 20, TRUE), labels=c("unverified", "verified")),
        wine = factor(sample(3, 20, TRUE), labels=c("red", "white", "none")),
        fruit = factor(sample(4, 20, TRUE), labels = c("apple", "banana", "cantaloup", "durian")),
        stars = factor(1:5, ordered=TRUE),
        tie_break = factor(1:5, 
                           labels=c("SetAlice", "AdvAlice", "Deuce", "AdvBob", "SetBob"),
                           ordered = TRUE)
    )
    y <- gl(2, 10)
    cv <- resample("crossvalidation", y)
    sets <- pre_split(df, y, cv[[1]])
    expect_is(pre_factor_to_logical(sets), "preprocessed_data")
    expect_is(pre_factor_to_logical(sets,
        base = c(wine = 3L, tie_break = 3L)))
    expect_is(pre_factor_to_logical(sets, 
        base = list(wine = 3L, fruit = NULL, tie_break = "Deuce")),
        "preprocessed_data")

    expect_error(pre_factor_to_logical(sets, feature="Dave Brubeck"))
    expect_error(pre_factor_to_logical(sets, base=c(wine=6L)))
    expect_error(pre_factor_to_logical(sets, base=c(fruit="xylophone")))
    expect_error(pre_factor_to_logical(sets, base=list(stars=NULL)))
    expect_warning(pre_factor_to_logical(sets, base=list(stars=1)))

    sets <- pre_factor_to_logical(sets,
                base = c(wine="none", tie_break="Deuce"),
                drop = c(fruit = FALSE, stars = FALSE, tie_break=FALSE))
    expect_identical(sets$feature_selection,
        structure(c(1L, 2L, 2L, 3L, 3L, 3L, 3L, 4L, 4L, 4L, 4L, 4L, 5L, 
        5L, 5L, 5L, 5L), .Names = c("email", "red", "white", "apple", 
        "banana", "cantaloup", "durian", "1", "2", "3", "4", "5", "SetAlice", 
        "AdvAlice", "Deuce", "AdvBob", "SetBob")))
})

