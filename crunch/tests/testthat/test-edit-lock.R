context("Locking and unlocking edit privileges")

if (run.integration.tests) {
    with(test.authentication, {
        with(test.dataset(df), {
            test_that("I can lock and unlock the dataset", {
                lock(ds)
                expect_error(name(ds) <- "Locked name",
                    "unlock")
                expect_error(categories(ds$v4) <- rev(categories(ds$v4)),
                    "unlock")
                expect_false(name(refresh(ds)) == "Locked name")
                unlock(ds)
                name(ds) <- "Unlocked name"
                expect_identical(name(refresh(ds)), "Unlocked name")
            })
        })
    })
}
