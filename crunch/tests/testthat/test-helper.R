context("Things in helper.R")

with(temp.option(foo.bar="no", foo.other="other"), {
    ## Mock outside of test_that because the expectations also check Sys.getenv
    with_mock(`Sys.getenv`=function (q, ...) {
            return(list(R_FOO_BAR="yes")[[q]] %||% "")
        },
        e1 <- envOrOption("foo.bar"),
        e2 <- envOrOption("foo.other"),
        e3 <- envOrOption("somethingelse")
    )

    test_that("envOrOption gets the right thing", {
        expect_identical(e1, "yes") ## Env var trumps option
        expect_identical(envOrOption("foo.bar"), "no") ## Bc no more env var
        expect_identical(e2, "other") ## Option if there is no env var
        expect_identical(e3, NULL) ## Null if neither
    })
})

test_that("SUTD", {
    a <- NULL
    tester <- setup.and.teardown(function () a <<- FALSE,
        function () a <<- TRUE)

    expect_true(is.null(a))
    with(tester, {
        expect_false(is.null(a))
        expect_false(a)
        ## Test that assertion failures are raised
        # expect_false(TRUE)
    })
    expect_true(a)

    ## Test that even if the code in the with block throws an error, (1) the
    ## teardown is run, and (2) it doesn't fail silently but turns into a
    ## failed test expectation.
    # a <- NULL
    # expect_true(is.null(a))
    # with(tester, {
    #     expect_false(is.null(a))
    #     expect_false(a)
    #     halt("Testing error handling, please ignore")
    # })
    # expect_true(a)
})

test_that("JSON equivalence", {
    expect_json_equivalent(4, 4)
    expect_json_equivalent(list(n=5, q=list(r=55, p=9)),
        list(q=list(p=9, r=55), n=5))
    expect_json_equivalent(list(1, 2), list(1, 2))
    skip("Upgrade these tests after testthat version lands")
    exp <- json_equivalent(list(1, 2))(list(2, 1))
    expect_false(exp$passed)
    expect_match(exp$failure, "not JSON-equivalent to")

    obj <- list(c=1, b=list(list(2, 3), list(d=9, f=5)), a=5)
    test <- json_equivalent(obj)
    expect_json_equivalent(list(b=list(list(2, 3), list(f=5, d=9)), c=1, a=5),
        obj)
    expect_false(test(list(c=1, b=list(list(3, 2), list(d=9, f=5)), a=5))$passed)
    expect_false(test(list(c=1, b=list(list(d=9, f=5), list(2, 3)), a=5))$passed)
})
