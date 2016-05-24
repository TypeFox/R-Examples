context("Shallow copies of variables")

with_mock_HTTP({
    ds <- loadDataset("test ds")
    expect_true(inherits(copy(ds$gender), "VariableDefinition"))
    expected <- VariableDefinition(
        name="Gender (copy)",
        alias="gender_copy",
        description="Gender",
        discarded=FALSE,
        format=list(summary=list(digits=2)),
        view=list(include_missing=FALSE,
            show_counts=FALSE,
            show_codes=FALSE,
            column_width=NULL
        ),
        expr=list(
            `function`="copy_variable",
            args=list(
                list(variable="/api/datasets/dataset1/variables/gender.json")
            )
        )
    )
    expect_json_equivalent(copy(ds$gender), expected)
})

if (run.integration.tests) {
    with(test.authentication, {
        with(test.dataset(newDatasetFromFixture("apidocs")), {
            q1_url <- self(ds$q1)
            varcat_url <- self(variables(ds))
            test_that("Can copy and manipulate a categorical variable", {
                expect_false("copy1" %in% names(ds))
                expect_true("q1" %in% names(ds))
                expect_silent(ds$copy1 <- copy(ds$q1, name="copy1"))
                expect_identical(as.vector(ds$copy1), as.vector(ds$q1))
                expect_false(name(ds$copy1) == name(ds$q1))
                expect_false(alias(ds$copy1) == alias(ds$q1))
                expect_false(self(ds$copy1) == self(ds$q1))
                ds <- refresh(ds)
                expect_true("copy1" %in% names(ds))
                expect_true("q1" %in% names(ds))

                skip("Can't edit categories in copy")
                ## Edit category in copy
                names(categories(ds$copy1))[2] <- "Canine"
                expect_identical(names(categories(ds$copy1))[1:3],
                    c("Cat", "Canine", "Bird"))
                expect_identical(names(categories(ds$q1))[1:3],
                    c("Cat", "Dog", "Bird"))

                ## Edit categories in original
                categories(ds$q1)[1:2] <- categories(ds$q1)[2:1]
                expect_identical(names(categories(ds$copy1))[1:3],
                    c("Cat", "Canine", "Bird"))
                expect_identical(names(categories(ds$q1))[1:3],
                    c("Dog", "Cat", "Bird"))
            })

            test_that("Can copy an array variable and manipulate it independently", {
                ds$allpets2 <- copy(ds$allpets)
                expect_true("allpets" %in% names(ds))
                expect_true("allpets2" %in% names(ds))
                expect_identical(name(ds$allpets2), "All pets owned (copy)")
                name(ds$allpets2) <- "Copy of allpets"
                expect_identical(name(ds$allpets2), "Copy of allpets")
                expect_identical(name(ds$allpets), "All pets owned")

                ## Edit subvariables in the copy
                subvariables(ds$allpets2)[1:2] <- subvariables(ds$allpets2)[2:1]
                expect_identical(names(subvariables(ds$allpets2)),
                    c("Dog", "Cat", "Bird"))
                expect_identical(names(subvariables(ds$allpets)),
                    c("Cat", "Dog", "Bird"))

                ## Edit subvariable names in the original
                names(subvariables(ds$allpets))[2] <- "Canine"
                expect_identical(names(subvariables(ds$allpets2)),
                    c("Dog", "Cat", "Bird"))
                expect_identical(names(subvariables(ds$allpets)),
                    c("Cat", "Canine", "Bird"))
            })

            test_that("Can copy subvariables (as non-subvars)", {

            })

            test_that("Can make a new array using copies of other array's subvars", {
                ## Let's make two copies of allpets, give them different names,
                ## then unbind them and rebind their subvars by animal
                ds$allpets_adult <- copy(ds$allpets, name="Adult pets")
                ds$allpets_juv <- copy(ds$allpets, name="Juvenile pets")
                unbind(ds$allpets_adult)
                unbind(ds$allpets_juv)
                ds <- refresh(ds)

                ds$allcats <- makeMR(pattern="Cat$", data=ds,
                    selections="selected", name="All cats", key="name")
                ds$alldogs <- makeMR(pattern="Canine$", data=ds,
                    selections="selected", name="All dogs", key="name")
                ds$allbirds <- makeMR(pattern="Bird$", data=ds,
                    selections="selected", name="All birds", key="name")

                expect_equivalent(as.array(crtabs(~ allcats, data=ds)),
                    array(c(4, 4), dim=c(2L),
                        dimnames=list(allcats=c("Adult pets | Cat",
                        "Juvenile pets | Cat"))))
            })

            with(test.dataset(newDatasetFromFixture("apidocs"), "part2"), {
                test_that("Copies get data when appending", {
                    out <- appendDataset(ds, part2)
                    expect_equivalent(as.array(crtabs(~ allcats, data=out)),
                        array(c(8, 8), dim=c(2L),
                            dimnames=list(allcats=c("Adult pets | Cat",
                            "Juvenile pets | Cat"))))
                })
            })
        })
    })
}
