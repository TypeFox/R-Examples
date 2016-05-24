context("Context manager")

test_that("Enter and exit are run", {
    a <- NULL
    tester <- ContextManager(function () a <<- FALSE,
        function () a <<- TRUE)

    expect_true(is.null(a))
    with(tester, {
        expect_false(is.null(a))
        expect_false(a)
        ## Test that assertion failures are raised
        # expect_false(TRUE)
    })
    expect_true(a)
})

test_that("An error in the code being executed is thrown but exit still runs", {
    a <- NULL
    tester <- ContextManager(function () a <<- FALSE,
        function () a <<- TRUE)

    expect_true(is.null(a))
    expect_error(with(tester, {
        halt("Testing error handling")
    }), "Testing error handling")
    expect_true(a)
})

test_that("Custom error handlers", {
    a <- NULL
    tester <- ContextManager(function () a <<- FALSE,
        function () a <<- TRUE,
        error=function (e) halt("Custom error"))

    expect_true(is.null(a))
    expect_error(with(tester, {
        halt("Testing error handling")
    }), "Custom error")
    expect_true(a)
})

test_that("'as' argument for output of enter function", {
    a <- FALSE
    ctx <- ContextManager(function () return(1:4),
        function () a <<- TRUE)

    expect_false(a)
    with(ctx, as="b", {
        expect_equivalent(sum(b), 10)
        d <- sum(b)
    })
    expect_true(a)
    expect_equivalent(d, 10)
})

test_that("'as' specified in the context manager itself", {
    a <- FALSE
    ctx <- ContextManager(function () return(1:4),
        function () a <<- TRUE,
        as="b")

    expect_false(a)
    with(ctx, {
        expect_equivalent(sum(b), 10)
    })
    expect_true(a)
})

test_that("temp.options", {
    options(crunch.test.option.test="foo")
    expect_identical(getOption("crunch.test.option.test"), "foo")
    expect_identical(getOption("crunch.test.test.test.test"), NULL)
    with(temp.options(crunch.test.option.test="bar",
                     crunch.test.test.test.test="test"), {
        expect_identical(getOption("crunch.test.option.test"), "bar")
        expect_identical(getOption("crunch.test.test.test.test"), "test")
    })
    expect_identical(getOption("crunch.test.option.test"), "foo")
    expect_identical(getOption("crunch.test.test.test.test"), NULL)
})

test_that("consent", {
    expect_false(askForPermission()) ## Because not interactive running
    with(consent(), {
        expect_true(askForPermission())
    })
    expect_false(askForPermission())
})
