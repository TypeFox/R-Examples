context("Categories")

with_mock_HTTP({
    ds <- loadDataset("test ds")
    cats <- categories(ds$gender)

    test_that("category init", {
        expect_true(inherits(cats[[1]], "Category"))
        expect_true(is.category(cats[[1]]))
        expect_true(inherits(cats, "Categories"))
        expect_true(is.categories(cats))
        expect_identical(length(cats), 3L)
    })

    test_that("Categories validation", {
        expect_error(Categories(
            list(id=-1L, name="B", numeric_value=1L, missing=FALSE),
            list(id=2L, name="C", numeric_value=2L, missing=FALSE),
            list(id=-1L, name="No Data", numeric_value=NULL, missing=TRUE)
        ), "Invalid category ids: must be unique")
        expect_error(Categories(
            list(id=1L, name="Name 1", numeric_value=1L, missing=FALSE),
            list(id=2L, name="Name 1", numeric_value=2L, missing=FALSE),
            list(id=-1L, name="No Data", numeric_value=NULL, missing=TRUE)
        ), "Invalid category names: must be unique")
    })

    test_that("category slicers", {
        expect_true(is.categories(cats[1]))
        expect_error(cats[c(1, 2, 5)],
            "subscript out of bounds: 5")
        expect_error(cats[c(1, 2, 98, 99)],
            "subscript out of bounds: 98 and 99")
    })

    test_that("can use negative subscripts on Categories", {
        expect_true(is.categories(cats[-1]))
        expect_error(cats[c(1, -1)],
            "only 0's may be mixed with negative subscripts")
    })

    test_that("categories to/fromJSON", {
        ## cereal serializes to JSON and then deserializes
        expect_identical(cats, Categories(data=cereal(cats)))
        expect_identical(cats[1], Categories(data=cereal(cats[1])))
        expect_identical(cats[[1]], Category(data=cereal(cats[[1]])))
    })

    test_that("lapply categories", {
        expect_identical(lapply(cats, function (x) x), cats)
    })

    test_that("category getters", {
        male <- cats[[1]]
        expect_identical(name(male), "Male")
        expect_identical(value(male), 1)
        expect_identical(id(male), 1L)
        expect_identical(names(cats), c("Male", "Female", "No Data"))
        expect_identical(values(cats), c(1, 2, NA))
        expect_identical(ids(cats), c(1L, 2L, -1L))
    })

    test_that("category setters", {
        male <- cats[[1]]
        name(male) <- "uomo"
        expect_identical(name(male), "uomo")
        expect_true(is.category(male))
        value(male) <- 42
        expect_identical(value(male), 42)
        expect_error(id(male) <- 4)
        expect_error(value(male) <- "foo")
        expect_true(is.category(male))
    })

    test_that("categories setters", {
        new_names <- c("masculino", "femenino", "No Data")
        names(cats) <- new_names
        expect_true(is.categories(cats))
        expect_equal(names(cats), new_names)
        names(cats)[2] <- "donne"
        expect_true(is.categories(cats))
        expect_equal(names(cats)[1:2], c("masculino", "donne"))
    })

    test_that("validation on category setting", {
        expect_error(cats[1] <- "new name",
            "Invalid categories: 1 element is not a Crunch category object")
        expect_error(name(cats[[1]]) <- NULL,
            'Names must be of class "character"')
    })

    test_that("names(categories)<- input validation", {
        expect_identical(names(categories(ds$gender)),
            c("Male", "Female", "No Data"))
        expect_error(names(categories(ds$gender))[2] <- NULL,
            "replacement has length zero") ## R default, good enough
        expect_error(names(categories(ds$gender))[2] <- list(foo="1"),
            'Names must be of class "character"')
        expect_error(names(categories(ds$gender))[23] <- "cat",
            "Invalid names: supplied 23 names for 3 categories")
            ## Not ideal error message, but best we can do here
    })

    test_that("categories ids cannot be set", {
        expect_error(ids(cats) <- rev(ids(cats)),
            "Cannot modify category ids")
    })

    test_that("dichotomize", {
        male <- cats[[1]]
        expect_false(is.selected(male))
        male$selected <- TRUE
        expect_true(is.selected(male))
        expect_equal(name(male), "Male")

        expect_false(is.dichotomized(cats))
        dcats <- dichotomize(cats, 1)
        expect_true(is.dichotomized(dcats))
        expect_true(is.selected(dcats[[1]]))
        expect_false(is.selected(dcats[[2]]))
        expect_equal(name(dcats[[1]]), "Male")
        expect_equal(names(dcats), c("Male", "Female", "No Data"))

        dcats2 <- dichotomize(cats, "Female")
        expect_true(is.dichotomized(dcats2))
        expect_false(is.selected(dcats2[[1]]))
        expect_true(is.selected(dcats2[[2]]))

        expect_error(dichotomize(cats, "Cat!"))

        cats2 <- undichotomize(dcats)
        expect_false(is.dichotomized(cats2))
        expect_false(is.selected(cats2[[1]]))
    })

    test_that("is.na", {
        expect_identical(is.na(cats), structure(c(FALSE, FALSE, TRUE),
            .Names=c("Male", "Female", "No Data")))
        expect_true(is.na(cats[[3]]))
        expect_false(is.na(cats[[1]]))
    })

    test_that("is.na<- by name", {
        cats <- cats
        try(is.na(cats) <- "Female")
        expect_true(is.categories(cats))
        expect_identical(is.na(cats), structure(c(FALSE, TRUE, TRUE),
            .Names=c("Male", "Female", "No Data")))
        expect_error(is.na(cats) <- c("Male", "Prefer not to say"),
            paste0("Category not found: ", dQuote("Prefer not to say")))
        expect_identical(is.na(cats), structure(c(FALSE, TRUE, TRUE),
            .Names=c("Male", "Female", "No Data")))
    })
    test_that("is.na<- by logical", {
        cats <- cats
        try(is.na(cats) <- c(TRUE, FALSE, FALSE))
        expect_true(is.categories(cats))
        expect_identical(is.na(cats), structure(c(TRUE, FALSE, FALSE),
            .Names=c("Male", "Female", "No Data")))
    })

    test_that("na.omit", {
        expect_identical(length(cats), 3L)
        expect_identical(length(na.omit(cats)), 2L)
        expect_true(is.categories(na.omit(cats)))
        expect_true(all(vapply(na.omit(cats), is.category, logical(1))))
    })

    newcat <- Category(name="Other", id=4)
    newcat2 <- Category(name="Something else", id=5)
    cats2 <- Categories(newcat, newcat2)
    test_that("Category constructor with missing attributes", {
        expect_false(is.na(newcat))
        expect_true(is.na(value(newcat)))
    })
    test_that("c() method for Categories, setup", {
        expect_true(is.categories(cats))
        expect_true(is.categories(cats2))
        expect_true(is.category(newcat))
    })
    test_that("c(Categories, Category)", {
        expect_true(is.categories(c(cats, newcat)))
    })
    test_that("c(Category, Categories)", {
        expect_true(is.categories(c(newcat, cats)))
    })
    test_that("c(Category, Category)", {
        expect_true(is.categories(c(newcat, newcat2)))
    })
    test_that("c(Categories, Categories)", {
        expect_true(is.categories(c(cats, cats2)))
    })
})



if (run.integration.tests) {
    with(test.authentication, {
        with(test.dataset(df), {
            test_that("categories setters persist to the server", {
                expect_equal(names(categories(ds$v4)), c("B", "C", "No Data"))
                names(categories(ds$v4))[1] <- "V"
                expect_equal(names(categories(ds$v4)), c("V", "C", "No Data"))
                expect_identical(names(categories(ds$v4)),
                    names(categories(refresh(ds)$v4)))

                categories(ds$v4)[1:2] <- categories(ds$v4)[2:1]
                expect_equal(names(categories(ds$v4)), c("C", "V", "No Data"))
            })

            test_that("categories<- with invalid input gives helpful message", {
                expect_error(categories(ds$v4) <- 1:3,
                    paste("`categories(x) <- value` only accepts Categories,",
                        "not numeric. Did you mean",
                        "`values(categories(x)) <- value`?"),
                    fixed=TRUE)
                expect_error(categories(ds$v4) <- c("A", "B", "C"),
                    paste("`categories(x) <- value` only accepts Categories,",
                        "not character. Did you mean",
                        "`names(categories(x)) <- value`?"),
                    fixed=TRUE)
                expect_error(categories(ds$v4) <- list(),
                    "`categories(x) <- value` only accepts Categories, not list.",
                    fixed=TRUE)
                expect_error(categories(ds$v1) <- 1:3,
                    "category assignment not defined for NumericVariable")
                expect_error(categories(ds$v4) <- categories(ds$v4)[c(1, 2, 5)],
                    "subscript out of bounds: 5")
            })
        })

        with(test.dataset(df), {
            test_that("Can add categories with c()", {
                expect_identical(names(categories(ds$v4)),
                    c("B", "C", "No Data"))
                categories(ds$v4) <- c(categories(ds$v4),
                    Category(name="D", id=4))
                expect_identical(names(categories(ds$v4)),
                    c("B", "C", "No Data", "D"))
            })
            test_that("Can insert a category in the middle", {
                ds$v4b <- df$v4
                expect_identical(names(categories(ds$v4b)),
                    c("B", "C", "No Data"))
                categories(ds$v4b) <- c(categories(ds$v4b)[1:2],
                    Category(name="D", id=4), categories(ds$v4b)[3])
                expect_identical(names(categories(ds$v4b)),
                    c("B", "C", "D", "No Data"))
            })
            test_that("Can add one to the end", {
                ds$v4c <- df$v4
                expect_identical(names(categories(ds$v4c)),
                    c("B", "C", "No Data"))
                categories(ds$v4c)[[4]] <- Category(name="D", id=4)
                expect_identical(names(categories(ds$v4c)),
                    c("B", "C", "No Data", "D"))
            })
            test_that("Can't duplicate categories", {
                ds$v4d <- df$v4
                expect_identical(names(categories(ds$v4d)),
                    c("B", "C", "No Data"))
                expect_error(categories(ds$v4) <- c(categories(ds$v4d),
                    categories(ds$v4d)))
            })
            test_that("Can delete a category that has no data", {
                ds$v4e <- df$v4
                categories(ds$v4e) <- c(categories(ds$v4e)[1:2],
                    Category(name="D", id=4), categories(ds$v4e)[3])
                expect_identical(names(categories(ds$v4e)),
                    c("B", "C", "D", "No Data"))
                ## Reassign the data from C to D
                ds$v4e[ds$v4e == "C"] <- "D"
                expect_equivalent(as.array(crtabs(~ v4e, data=ds)),
                    array(c(10, 0, 10), dim=3L, dimnames=list(v4=c("B", "C", "D"))))
                ## Then delete B
                categories(ds$v4e) <- categories(ds$v4e)[-2]
                expect_identical(names(categories(ds$v4e)),
                    c("B", "D", "No Data"))
                expect_equivalent(as.array(crtabs(~ v4e, data=ds)),
                    array(c(10, 10), dim=2L, dimnames=list(v4=c("B", "D"))))
            })
        })

        with(test.dataset(mrdf), {
            ds <- mrdf.setup(ds)
            test_that("Can edit and reorder categories in categorical array", {
                expect_identical(names(categories(ds$CA)),
                    c("0.0", "1.0", "No Data"))
                names(categories(ds$CA))[1:2] <- c("First", "Second")
                expect_identical(names(categories(ds$CA)),
                    c("First", "Second", "No Data"))
                categories(ds$CA) <- rev(categories(ds$CA))
                expect_identical(names(categories(ds$CA)),
                    c("No Data", "Second", "First"))
            })
        })

        with(test.dataset(newDatasetFromFixture("apidocs")), {
            test_that("dichotomizing dichotomizes the subvariables", {
                expect_true(is.MR(ds$allpets))
                expect_true(is.dichotomized(categories(ds$allpets)))
                expect_true(is.dichotomized(categories(ds$allpets$Dog)))
                ds$allpets <- undichotomize(ds$allpets)
                expect_true(is.CA(ds$allpets))
                expect_false(is.dichotomized(categories(ds$allpets)))
                expect_false(is.dichotomized(categories(ds$allpets$Dog)))
            })
            test_that("Editing array categories affects the subvariables too", {
                expect_identical(names(categories(ds$petloc)),
                    c("Cat", "Dog", "Bird", "Skipped", "Not Asked"))
                expect_identical(names(categories(ds$petloc$Home)),
                    c("Cat", "Dog", "Bird", "Skipped", "Not Asked"))
                expect_identical(names(table(ds$petloc$Home)),
                    c("Cat", "Dog", "Bird"))
                names(categories(ds$petloc))[2] <- "Canine"
                expect_identical(names(categories(ds$petloc)),
                    c("Cat", "Canine", "Bird", "Skipped", "Not Asked"))
                expect_identical(names(table(ds$petloc$Home)),
                    c("Cat", "Canine", "Bird"))
                expect_identical(names(categories(ds$petloc$Home)),
                    c("Cat", "Canine", "Bird", "Skipped", "Not Asked"))
            })
        })

        with(test.dataset(df), {
            test_that("Cache invalidation when modifying categories", {
                expect_equal(names(categories(ds$v4)), c("B", "C", "No Data"))
                expect_equivalent(as.array(crtabs(~ v4, data=ds)),
                    array(c(10, 10), dim=2L, dimnames=list(v4=c("B", "C"))))

                names(categories(ds$v4))[1] <- "V"
                expect_equal(names(categories(ds$v4)), c("V", "C", "No Data"))
                expect_equivalent(as.array(crtabs(~ v4, data=ds)),
                    array(c(10, 10), dim=2L, dimnames=list(v4=c("V", "C"))))
            })
        })
    })
}
