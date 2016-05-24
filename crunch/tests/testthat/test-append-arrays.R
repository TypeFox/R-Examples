context("Appending datasets with arrays")

if (run.integration.tests) {
    with(test.authentication, {
        with(test.dataset(mrdf, "part1"), {
            part1 <- mrdf.setup(part1, selections="1.0")
            with(test.dataset(mrdf, "part2"), {
                part2 <- mrdf.setup(part2, selections="1.0")
                test_that("set up MR for appending", {
                    expect_true(is.Multiple(part1$MR))
                    expect_equivalent(as.array(crtabs(~ MR, data=part1)),
                        array(c(2, 1, 1), dim=c(3L),
                        dimnames=list(MR=c("mr_1", "mr_2", "mr_3"))))
                    expect_true(is.Multiple(part2$MR))
                    expect_equivalent(as.array(crtabs(~ MR, data=part2)),
                        array(c(2, 1, 1), dim=c(3L),
                        dimnames=list(MR=c("mr_1", "mr_2", "mr_3"))))
                    # 2 because there is always a "batch 0" which
                    # is present for rows which were added directly
                    # rather than through a batch append,
                    # and another for "part1".
                    expect_identical(length(batches(part1)), 2L)
                    expect_identical(length(batches(part2)), 2L)
                })
                test_that("identical datasets with arrays can append", {
                    expect_message(out <- appendDataset(part1, part2),
                        "No conflicts")
                    expect_true(is.dataset(out))
                    expect_identical(length(batches(out)), 3L)
                    expect_identical(dim(out), c(nrow(mrdf)*2L, 2L))
                    expect_true(is.Multiple(out$MR))
                    expect_equivalent(as.array(crtabs(~ MR, data=out)),
                        array(c(4, 2, 2), dim=c(3L),
                        dimnames=list(MR=c("mr_1", "mr_2", "mr_3"))))
                    expect_true(is.dataset(refresh(part2)))
                })
            })
        })

        with(test.dataset(mrdf, "part1"), {
            part1 <- mrdf.setup(part1, selections="1.0")
            mr_cats <- categories(part1$MR)
            subvar_cats <- categories(part1$MR$mr_1)
            dichotomized_cats <- Categories(
                list(id=2L, missing=FALSE, name="0.0", numeric_value=0, selected=FALSE),
                list(id=1L, missing=FALSE, name="1.0", numeric_value=1, selected=TRUE),
                list(id=-1L, missing=TRUE, name="No Data", numeric_value=NULL, selected=FALSE))
            with(test.dataset(mrdf, "part2"), {
                ## Dichotomize this way so that categories get aligned
                ## (via supertype)
                part2 <- mrdf.setup(part2)
                unbind(part2$CA)
                part2 <- refresh(part2)
                undichotomized_cats <- Categories(
                    list(id=2L, missing=FALSE, name="0.0", numeric_value=0),
                    list(id=1L, missing=FALSE, name="1.0", numeric_value=1),
                    list(id=-1L, missing=TRUE, name="No Data", numeric_value=NULL))
                test_that("set up MR for appending", {
                    expect_true(is.Multiple(part1$MR))
                    expect_equivalent(as.array(crtabs(~ MR, data=part1)),
                        array(c(2, 1, 1), dim=c(3L),
                        dimnames=list(MR=c("mr_1", "mr_2", "mr_3"))))
                    expect_true(is.null(part2$MR))
                    expect_identical(mr_cats, subvar_cats)
                    expect_identical(mr_cats, dichotomized_cats)
                    expect_identical(categories(part2$mr_1),
                        undichotomized_cats)
                    expect_false(identical(dichotomized_cats,
                        undichotomized_cats)) ## Just being clear about that
                    expect_identical(as.vector(part1$MR$mr_1),
                        as.vector(part2$mr_1))
                    expect_identical(as.vector(part1$MR$mr_2),
                        as.vector(part2$mr_2))
                    expect_identical(as.vector(part1$MR$mr_3),
                        as.vector(part2$mr_3))
                })
                out <- suppressMessages(try(appendDataset(part1, part2)))
                test_that("Dataset #2 isn't modified by appending to another", {
                    expect_identical(refresh(part2), part2)
                })
                test_that("unbound subvariables get lined up", {
                    expect_true(is.dataset(out))
                    expect_identical(length(batches(out)), 3L)
                    expect_identical(dim(out), c(nrow(mrdf)*2L, 2L))
                    expect_true(is.variable(out$MR))
                    expect_identical(categories(out$MR), dichotomized_cats)
                    expect_identical(categories(out$MR$mr_1), dichotomized_cats)
                    expect_false(identical(categories(out$MR),
                        undichotomized_cats))
                    expect_identical(as.vector(out$MR$mr_1),
                        rep(as.vector(part2$mr_1), 2))
                    expect_true(is.Multiple(out$MR))
                    expect_identical(names(subvariables(out$MR)),
                        c("mr_1", "mr_2", "mr_3"))
                    expect_equivalent(as.array(crtabs(~ MR, data=out)),
                        array(c(4, 2, 2), dim=c(3L),
                        dimnames=list(MR=c("mr_1", "mr_2", "mr_3"))))
                })
            })
        })
        with(test.dataset(mrdf, "part1"), {
            part1 <- mrdf.setup(part1, selections="1.0")
            mr_cats <- categories(part1$MR)
            subvar_cats <- categories(part1$MR$mr_1)
            dichotomized_cats <- Categories(
                list(id=2L, missing=FALSE, name="0.0", numeric_value=0, selected=FALSE),
                list(id=1L, missing=FALSE, name="1.0", numeric_value=1, selected=TRUE),
                list(id=-1L, missing=TRUE, name="No Data", numeric_value=NULL, selected=FALSE))
            with(test.dataset(mrdf, "part2"), {
                cast.these <- grep("mr_", names(part2))
                part2[cast.these] <- lapply(part2[cast.these],
                    castVariable, "categorical")
                undichotomized_cats <- Categories(
                    list(id=2L, missing=FALSE, name="0.0", numeric_value=0),
                    list(id=1L, missing=FALSE, name="1.0", numeric_value=1),
                    list(id=-1L, missing=TRUE, name="No Data", numeric_value=NULL))
                test_that("set up MR for appending", {
                    expect_true(is.Multiple(part1$MR))
                    expect_equivalent(as.array(crtabs(~ MR, data=part1)),
                        array(c(2, 1, 1), dim=c(3L),
                        dimnames=list(MR=c("mr_1", "mr_2", "mr_3"))))
                    expect_true(is.null(part2$MR))
                    expect_identical(mr_cats, subvar_cats)
                    expect_identical(mr_cats, dichotomized_cats)
                    expect_identical(categories(part2$mr_1),
                        undichotomized_cats)
                    expect_identical(as.vector(part1$MR$mr_1),
                        as.vector(part2$mr_1))
                    expect_identical(as.vector(part1$MR$mr_2),
                        as.vector(part2$mr_2))
                    expect_identical(as.vector(part1$MR$mr_3),
                        as.vector(part2$mr_3))
                })
                test_that("unbound subvars with not identical cats", {
                    expect_message(out <- appendDataset(part1, part2),
                        "No conflicts")
                    expect_true(is.dataset(out))
                    expect_identical(length(batches(out)), 3L)
                    expect_identical(dim(out), c(nrow(mrdf)*2L, 2L))
                    expect_true(is.variable(out$MR))
                    expect_identical(categories(out$MR), dichotomized_cats)
                    expect_identical(categories(out$MR$mr_1), dichotomized_cats)
                    expect_false(identical(categories(out$MR),
                        undichotomized_cats)) ## To be clear about the problem
                    ## Not handling categories with different ids but same names
                    expect_identical(as.vector(out$MR$mr_1),
                        rep(as.vector(part2$mr_1), 2))
                    expect_true(is.Multiple(out$MR))
                    expect_identical(names(subvariables(out$MR)),
                        c("mr_1", "mr_2", "mr_3"))
                    expect_equivalent(as.array(crtabs(~ MR, data=out)),
                        array(c(4, 2, 2), dim=c(3L),
                        dimnames=list(MR=c("mr_1", "mr_2", "mr_3"))))
                })
            })
        })

        with(test.dataset(mrdf[-3], "part1"), {
            part1 <- mrdf.setup(part1, selections="1.0")
            part1 <- saveVersion(part1, "Before appending")
            with(test.dataset(mrdf[-1], "part2"), {
                part2 <- mrdf.setup(part2, selections="1.0")
                test_that("set up MR for appending", {
                    expect_true(is.Multiple(part1$MR))
                    expect_identical(names(subvariables(part1$MR)),
                        c("mr_1", "mr_2"))
                    expect_equivalent(as.array(crtabs(~ MR, data=part1)),
                        array(c(2, 1), dim=c(2L),
                        dimnames=list(MR=c("mr_1", "mr_2"))))
                    expect_true(is.Multiple(part2$MR))
                    expect_identical(names(subvariables(part2$MR)),
                        c("mr_2", "mr_3"))
                    expect_equivalent(as.array(crtabs(~ MR, data=part2)),
                        array(c(1, 1), dim=c(2L),
                        dimnames=list(MR=c("mr_2", "mr_3"))))
                })
                out <- suppressMessages(try(appendDataset(part1, part2)))
                test_that("arrays with different subvariables can append", {
                    expect_true(is.dataset(out))
                    expect_identical(length(batches(out)), 3L)
                    expect_identical(dim(out), c(nrow(mrdf)*2L, 2L))
                    expect_true(is.variable(out$MR))
                    expect_true(is.Multiple(out$MR))
                    expect_identical(names(subvariables(out$MR)),
                        c("mr_1", "mr_2", "mr_3"))
                    expect_equivalent(as.array(crtabs(~ MR, data=out)),
                        array(c(2, 2, 1), dim=c(3L),
                        dimnames=list(MR=c("mr_1", "mr_2", "mr_3"))))
                })

                test_that("Rolling back to initial import reverts the append", {
                    out <- restoreVersion(out, "Before appending")
                    expect_true(is.Multiple(out$MR))
                    expect_identical(names(subvariables(out$MR)),
                        c("mr_1", "mr_2"))
                    expect_equivalent(as.array(crtabs(~ MR, data=out)),
                        array(c(2, 1), dim=c(2L),
                        dimnames=list(MR=c("mr_1", "mr_2"))))
                    expect_identical(length(batches(out)), 2L)
                })
            })
        })

        sparsemrdf1 <- data.frame(v4=factor(rep(c("a", "b"), 500)))
        sparsemrdf2 <- data.frame(mr_1=c(1,0,1,1,0, rep(NA, 995)),
                           mr_2=c(rep(NA, 995), 0, 1, 1, 1, 0),
                           v4=as.factor(LETTERS[2:3]))
        with(test.dataset(sparsemrdf1, "part1"), {
            with(test.dataset(sparsemrdf2, "part2"), {
                part2 <- mrdf.setup(part2)
                out <- suppressMessages(try(appendDataset(part1, part2)))
                test_that("Sparse append with array", {
                    expect_identical(length(batches(out)), 3L)
                    expect_identical(nrow(out), 2000L)
                    expect_identical(as.vector(out$CA$mr_2),
                        factor(c(rep(NA, 1995), "0.0", "1.0", "1.0", "1.0", "0.0")))
                })

                test_that("Rolling back to initial import reverts the append", {
                    out <- restoreVersion(out, length(versions(out))) ## Get the oldest
                    expect_identical(nrow(out), 1000L)
                    expect_identical(length(batches(out)), 2L)
                })
            })
        })

        with(test.dataset(mrdf, "part1"), {
            part1 <- mrdf.setup(part1, selections="1.0")
            names(subvariables(part1$MR)) <- c("One", "Two", "Three")
            ## Aliases are c("mr_1", "mr_2", "mr_3")
            with(test.dataset(mrdf, "part2"), {
                part2 <- mrdf.setup(part2, selections="1.0")
                names(subvariables(part2$MR)) <- c("Loneliest", "Two", "Three")
                aliases(subvariables(part2$MR))[3] <- "alt"
                ## Aliases are c("mr_1", "mr_2", "alt")
                test_that("alias and name matching on appending arrays", {
                    expect_message(out <- appendDataset(part1, part2),
                        "No conflicts")
                    expect_true(is.dataset(out))
                    expect_identical(length(batches(out)), 3L)
                    expect_identical(dim(out), c(nrow(mrdf)*2L, 2L))
                    expect_true(is.Multiple(out$MR))
                    skip("We get 2 'Threes'")
                    expect_equivalent(as.array(crtabs(~ MR, data=out)),
                        array(c(4, 2, 2), dim=c(3L),
                        dimnames=list(MR=c("One", "Two", "Three"))))
                })
            })
        })

        with(test.dataset(newDatasetFromFixture("apidocs")), as="part1", {
            test_that("Setup for testing references post append", {
                expect_true(name(part1$allpets) == "All pets owned")
                name(part1$allpets) <- "Some of my pets"
                expect_true(name(part1$allpets) == "Some of my pets")
            })

            ## Release and re-lease
            part1 <- releaseAndReload(part1)

            test_that("Check again", {
                expect_true(name(part1$allpets) == "Some of my pets")
            })

            with(test.dataset(newDatasetFromFixture("apidocs")), as="part2", {
                out <- suppressMessages(try(appendDataset(part1, part2)))
                test_that("Append doesn't revert metadata changes", {
                    expect_false(name(out$allpets) == "All pets owned")
                    expect_true(name(out$allpets) == "Some of my pets")
                })

                ## Release and re-lease
                out <- releaseAndReload(out)

                test_that("Metadata sticks after releasing", {
                    expect_false(name(out$allpets) == "All pets owned")
                    expect_true(name(out$allpets) == "Some of my pets")
                })

                ## Change the name and release again
                name(out$allpets) <- "Apple"
                out <- releaseAndReload(out)
                test_that("Metadata sticks after releasing and not appending", {
                    expect_true(name(out$allpets) == "Apple")
                })
            })
        })
    })
}
