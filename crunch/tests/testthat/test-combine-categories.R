context("Combine categories")

with_mock_HTTP({
    ds <- loadDataset("test ds")

    both <- VariableDefinition(
        name="Gender 1 cat",
        description="Gender",
        discarded=FALSE,
        format=list(summary=list(digits=2)),
        view=list(include_missing=FALSE,
            show_counts=FALSE,
            show_codes=FALSE,
            column_width=NULL
        ),
        expr=list(
            `function`="combine_categories",
            args=list(
                list(variable="/api/datasets/dataset1/variables/gender.json"),
                list(value=list(
                    list(
                        name="Both",
                        combined_ids=I(c(1, 2)),
                        missing=FALSE,
                        numeric_value=NULL,
                        id=1
                    ),
                    list(
                        numeric_value=NULL,
                        missing=TRUE,
                        id=-1,
                        name="No Data",
                        combined_ids=I(-1)
                    )
                ))
            )
        )
    )
    test_that("combine() constructs the correct VarDef for categorical", {
        combine.names <- combine(ds$gender, name="Gender 1 cat",
            list(list(name="Both", categories=c("Male", "Female"))))
        expect_json_equivalent(combine.names, both)
        expect_json_equivalent(combine.names$expr, both$expr)
        combine.ids <- combine(ds$gender, name="Gender 1 cat",
            list(list(name="Both", categories=c(1,2))))
        expect_json_equivalent(combine.ids, both)
    })

    test_that("Default variable name for combine()", {
        expect_identical(combine(ds$gender,
            list(list(name="Both", categories=c("Male", "Female"))))$name,
            "Gender (1 category)")
        expect_identical(combine(ds$gender)$name,
            "Gender (2 categories)")
    })

    test_that("combine() validation on variable type", {
        expect_error(combine(ds$birthyr),
            paste0("Cannot combine ", dQuote("Birth Year"), ": must be type ",
            "categorical, categorical_array, or multiple_response"))
        expect_error(combine(ds$starttime),
            "categorical, categorical_array, or multiple_response")
    })

    test_that("combine() requires combinations", {
        expect_error(combine(ds$gender, combinations=45),
            "'combinations' must be a list of combination specifications")
        expect_error(combine(ds$gender,
            combinations=list(name="Both", categories=c("Male", "Female"))),
            "'combinations' must be a list of combination specifications")
        expect_error(combine(ds$gender,
            combinations=list(list(name="Both"))),
            "'combinations' must be a list of combination specifications")
    })

    test_that("combinations must have valid names", {
        expect_error(combine(ds$gender,
            list(list(name="Man", categories=1),
                list(name="Man", categories=2))),
            paste("Duplicate category name given:", dQuote("Man")))
        expect_error(combine(ds$gender,
            list(list(name="Male", categories=2))),
            paste("Duplicate category name given:", dQuote("Male")))
        expect_error(combine(ds$gender,
            list(list(name="Male", categories=1))),
            NA) ## "Male" is category in original but not result
    })

    test_that("combinations reference unique categories", {
        expect_error(combine(ds$gender, name="Gender 1 cat",
            list(
                list(name="Both", categories=c("Male", "Female")),
                list(name="Men", categories=1)
            )),
            paste("Category", dQuote("Male"),
                "referenced in multiple combinations"))
        expect_error(combine(ds$gender, name="Gender 1 cat",
            list(
                list(name="Both", categories=c("Male", "Female")),
                list(name="Men", categories=1),
                list(name="Ladies", categories="Female")
            )),
            paste("Categories", dQuote("Male"), "and", dQuote("Female"),
                "referenced in multiple combinations"))
    })

    test_that("combinations reference actual categories", {
        expect_error(combine(ds$gender, name="Gender 1 cat",
            list(list(name="Both", categories=c(1, 42)))),
            paste("Combination", dQuote("Both"),
                "references category with id 42, which does not exist"))
        expect_error(combine(ds$gender, name="Gender 1 cat",
            list(list(name="Both", categories=c("Male", "Not male")))),
            paste("Category not found:", dQuote("Not male")))
    })

    test_that("combinations reference categories by name or id (char or numeric)", {
        expect_error(combine(ds$gender, name="Gender 1 cat",
            list(list(name="Both", categories=list("Male", "Female")))),
            "Combinations must reference 'categories' by name or id")
        expect_error(combine(ds$gender, name="Gender 1 cat",
            list(list(name="Both", categories=NULL))),
            "Combinations must reference 'categories' by name or id")
    })
})

if (run.integration.tests) {
    with(test.authentication, {
        with(test.dataset(newDatasetFromFixture("apidocs")), {
            test_that("We can create a new categorical by combining", {
                ds$combined_pets <- combine(ds$q1, name="Pets (combined)",
                    list(list(name="Mammals", categories=c("Cat", "Dog"))))
                expect_identical(names(categories(ds$combined_pets)),
                    c("Mammals", "Bird", "Skipped", "Not Asked"))
                expect_equivalent(as.array(crtabs(~ combined_pets, data=ds)),
                    array(c(10, 3), dim=2,
                    dimnames=list(combined_pets=c("Mammals", "Bird"))))
            })

            test_that("combine() with no combinations is effectively a copy", {
                ds$combined_pets2 <- combine(ds$q1)
                expect_identical(as.vector(ds$combined_pets2), as.vector(ds$q1))
            })

            test_that("combine() with categorical array", {
                ds$combined_petloc <- combine(ds$petloc,
                    name="Pet locations (combined)",
                    list(list(name="Mammals", categories=c("Cat", "Dog"))))
                expect_identical(names(categories(ds$combined_petloc)),
                    c("Mammals", "Bird", "Skipped", "Not Asked"))
            })
        })
    })
}
