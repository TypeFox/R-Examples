context("Various helper functions")

test_that("is.error", {
    e <- try(halt("error in a box"), silent=TRUE)
    expect_true(is.error(e))
    expect_false(is.error("not an error"))
    expect_false(is.error(NULL))
    expect_false(is.error(NA))
    expect_error("not an error", NA)
})

test_that("update list", {
    a <- list(a=1, b=2)
    b <- list(c=3, b=4)
    expect_identical(updateList(a, b), list(a=1, b=4, c=3))
    expect_identical(updateList(list(), b), b)
    expect_identical(updateList(NULL, b), b)
})

test_that("selectFrom selects what it should", {
    l1 <- list(list(a=1, b=2), list(c=3, b=4))
    expect_identical(selectFrom("b", l1), c(2, 4))
    expect_identical(selectFrom("a", l1), c(1, NA))
    expect_identical(selectFrom("a", l1, ifnot=4), c(1, 4))
    expect_identical(selectFrom("d", l1), c(NA, NA))
    l2 <- l1
    l2[[2]] <- 4
    expect_identical(selectFrom("b", l2), c(2, NA))
    expect_error(selectFrom("b", 5), "xlist must be a list object")
})

test_that("rethrow a caught error", {
    e <- try(halt("error in a box"), silent=TRUE)
    expect_true(is.error(e))
    expect_error(rethrow(e), "error in a box")
})

test_that("%||%", {
    expect_identical("f" %||% "g", "f")
    expect_identical(NULL %||% "g", "g")
    expect_identical("f" %||% halt("Nooooooo!"), "f")
})

test_that("dirtyElements", {
    x <- list(
        list(a=1, b=1),
        list(a="1", b="1"),
        list(a="d", b="e")
    )
    y <- x
    expect_false(any(dirtyElements(x, y)))
    y[[2]]$b <- "f"
    y[[1]]$b <- 1
    expect_identical(dirtyElements(x, y), c(FALSE, TRUE, FALSE))
    y[[3]]$a <- "f"
    expect_identical(dirtyElements(x, y), c(FALSE, TRUE, TRUE))
})


test_that("joinPath", {
    expect_identical(joinPath("/api/datasets/", "../variables/"),
        "/api/variables/")
    expect_identical(joinPath("/api/variables/", "4412es.json"),
        "/api/variables/4412es.json")
    expect_identical(joinPath("a/b/c/d/../e/f/", "g/../../h/"),
        "a/b/c/e/h/")
    expect_identical(joinPath("/api/datasets/", "/variables/"),
        "/variables/")
    expect_identical(joinPath("/api/datasets/", "/"),
        "/")
    expect_identical(joinPath("/api/datasets/", "./id/"),
        "/api/datasets/id/")
})

test_that("absoluteURL", {
    base.url <- "https://fake.crunch.io/api/datasets/"
    expect_identical(absoluteURL("../variables/", base.url),
        "https://fake.crunch.io/api/variables/")
    expect_identical(absoluteURL("4412es.json", base.url),
        "https://fake.crunch.io/api/datasets/4412es.json")
    expect_identical(absoluteURL("g/../../h/",
        "https://fake.crunch.io/a/b/c/d/../e/f/"),
        "https://fake.crunch.io/a/b/c/e/h/")
    expect_identical(absoluteURL("/variables/", base.url),
        "https://fake.crunch.io/variables/")
    expect_identical(absoluteURL("/", base.url),
        "https://fake.crunch.io/")
})

test_that("emptyObject JSONifies correctly", {
    expect_equivalent(unclass(toJSON(emptyObject())), "{}")
})

test_that("setIfNotAlready", {
    with(temp.options(crunch.test.opt1="previous",
                      crunch.test.opt2=NULL,
                      crunch.test.opt3=4), {

        old <- setIfNotAlready(crunch.test.opt1="value", crunch.test.opt2=5)
        expect_identical(getOption("crunch.test.opt1"), "previous")
        expect_identical(getOption("crunch.test.opt2"), 5)
        expect_identical(getOption("crunch.test.opt3"), 4)
    })
})
