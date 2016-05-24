context("Debugging")

test_that("Debugging flag persistance", {
    if(R.Version()$major < 3 || R.Version()$minor < 1.0){
        is_debugged <- function(x){
            expectation(isdebugged(x), "isn't being debugged",
                        "is being debugged")
        }

        p <- modeling_procedure("lda")
        debug(p$fit_fun)
        expect_that(p$fit_fun, is_debugged)

        p <- listify(p)
        expect_that(p[[1]]$fit_fun, is_debugged)

        flags <- get_debug_flags(p)
        p[[1]]$importance_fun <- "Long ago before the world was cool the moon flew out of the earth."
        expect_false(
            isdebugged(p[[1]]$fit_fun)
        )
        expect_that(set_debug_flags(p, flags)[[1]]$fit_fun, is_debugged)
    }
})
