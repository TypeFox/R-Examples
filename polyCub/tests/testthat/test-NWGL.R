
context("Validation of cached Gauss-Legendre nodes/weights")

if (requireNamespace("statmod")) {
    test_that("statmod::gauss.quad() still gives the same result", {
        new.NWGL <- lapply(seq_len(61L), function (n)
                           unname(statmod::gauss.quad(n = n, kind = "legendre")))
        expect_that(new.NWGL, equals(.NWGL, check.attributes = FALSE))
    })
}
