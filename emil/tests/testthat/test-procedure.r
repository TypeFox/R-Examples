context("Modeling procedure")

test_that("Procedure coercing", {
    my_test <- function(x){
        expect_is(x, "list")
        expect_false(inherits(x, "modeling_procedure"))
        expect_true(all(sapply(x, inherits, "modeling_procedure")))
        expect_false(is.null(names(x)))
        expect_false(any(names(x) %in% c("", NA)))
    }
    my_test(multify("lda"))
    my_test(multify(c("qda", "lda")))
    my_test(multify(list("qda", "lda")))
    proc <- multify(c(a = "qda", "lda"))
    my_test(proc)
    expect_equal(names(proc), c("a", "lda"))
    my_test(multify(modeling_procedure("lda")))
    my_test(multify(list(modeling_procedure("lda"),
                         "qda")))

    expect_false(attr("lda", "multiple"))
    expect_true(attr(c("lda", "qda"), "multiple"))
    expect_true(attr(list("lda"), "multiple"))
})
