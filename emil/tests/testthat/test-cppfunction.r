context("Rcpp function")

test_that("Detection of constants", {
    f <- function(x, type, na.rm) is_constant(as(x, type), na.rm)
    for(p in 1:4){
        x <- as.matrix(do.call(expand.grid, rep(list(c(NA, 0, 1)), p)))
        is.constant <- c(NA, TRUE, FALSE)[sapply(apply(x, 1, table), length) + 1]
        is.constant.without.NA <- ifelse(is.constant,
            ifelse(apply(is.na(x), 1, any), NA, TRUE), FALSE)
        for(my_type in c("character", "complex", "expression", "factor", "integer", "logical", "numeric")){
            expect_identical(apply(x, 1, f, my_type, TRUE), is.constant)
            expect_identical(apply(x, 1, f, my_type, FALSE), is.constant.without.NA)
        }
    }
})

