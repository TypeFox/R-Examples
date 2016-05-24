context("Displaying append conflicts")

c1 <- list()
c2 <- list(
    var1=list(
        conflicts=list(list(
            message="No good",
            resolution="But I fixed it already"
        )),
        metadata=list(
            references=list(
                name="First"
            )
        )
    ),
    var4=list(conflicts=list())
)
c3 <- list(
    var2=list(
        conflicts=list(list(
            message="No good",
            resolution="But I fixed it already"
        )),
        metadata=list(
            references=list(
                name="Second"
            )
        )
    ),
    var4=list(conflicts=list()),
    var1=list(
        conflicts=list(
            list(
                message="No good",
                resolution="But I fixed it already"
            ),
            list(
                message="Oh, and there was another problem",
                resolution="But it's also cool"
            )
        ),
        metadata=list(
            references=list(
                name="First"
            )
        )
    )
)
c4 <- c3
c4$var3 <- list(
    conflicts=list(list(
        message="Type mismatch"
    )),
    metadata=list(references=list(name="Last"))
)

test_that("Simple conflict messages are formatted correctly", {
    expect_equivalent(flattenConflicts(c3),
        data.frame(
            message=c("No good", "No good", "Oh, and there was another problem"),
            resolution=c("But I fixed it already", "But I fixed it already", "But it's also cool"),
            url=c("var2", "var1", "var1"),
            name=c("Second", "First", "First"),
            stringsAsFactors=FALSE))
        expect_equivalent(flattenConflicts(c4),
            data.frame(
                message=c("No good", "No good", "Oh, and there was another problem", "Type mismatch"),
                resolution=c("But I fixed it already", "But I fixed it already", "But it's also cool", NA),
                url=c("var2", "var1", "var1", "var3"),
                name=c("Second", "First", "First", "Last"),
                stringsAsFactors=FALSE))

    expect_identical(formatConflicts(flattenConflicts(c1)), "No conflicts.")
    expect_identical(formatConflicts(flattenConflicts(c2)),
        paste("Conflict: No good; Resolution: But I fixed it already; 1 variable:", dQuote("First")))
    expect_identical(formatConflicts(flattenConflicts(c3)),
        c(paste("Conflict: No good; Resolution: But I fixed it already; 2 variables:", dQuote("Second"), "and", dQuote("First")),
        paste("Conflict: Oh, and there was another problem; Resolution: But it's also cool; 1 variable:", dQuote("First"))))
    expect_identical(formatFailures(flattenConflicts(c4)),
        paste("Critical conflict: Type mismatch; 1 variable:", dQuote("Last")))
})

source("conflicts.R")
test_that("Complex conflicts are formatted", {
    expect_identical(formatConflicts(flattenConflicts(mock.conflicts)),
        c(paste("Conflict: Only in existing dataset; Resolution: Additional rows will be marked missing.; 1 variable:", dQuote("mr_1")),
        paste("Conflict: Only in new dataset; Resolution: Variable will be added with existing rows marked missing.; 1 variable:", dQuote("mr_3")),
        paste("Conflict: Subvariables didn't match; Resolution: Union of subvariables will be used; 1 variable:", dQuote("MR"))))
})

if (run.integration.tests) {
    with(test.authentication, {
        with(test.dataset(mrdf, "part1"), {
            part1 <- mrdf.setup(part1)
            with(test.dataset(mrdf[c("mr_3", "v4")], "part2"), {
                alias(part2$mr_3) <- "CA"
                name(part2$CA) <- "Bad var"
                test_that("setup for append type mismatch", {
                    p1.batches <- batches(part1)
                    expect_true(inherits(p1.batches, "ShojiCatalog"))
                    expect_identical(length(p1.batches), 2L)
                    expect_true("CA" %in% names(part1))
                    expect_true("CA" %in% names(part2))
                    expect_true(is.CA(part1$CA))
                    expect_true(is.Numeric(part2$CA))
                })
                test_that("append conflict on type mismatch", {
                    expect_message(try(appendDataset(part1, part2),
                        silent=TRUE),
                        paste("Critical conflict: Variable is not array variable on both frames;",
                        "1 variable:", dQuote("Bad var")))
                    expect_error(appendDataset(part1, part2),
                        "There are conflicts that cannot be resolved automatically.")
                    expect_identical(length(batches(part1)), 2L)                })
            })
        })


    })
}
