context("User stuff")

with_mock_HTTP({
    test_that("Getting user object", {
        user <- getUser("/api/users/user1.json")
        expect_true(inherits(user, "ShojiObject"))
        expect_identical(user@body$email, "fake.user@example.com")
    })

    test_that("Getting account's user catalog", {
        usercat <- getAccountUserCatalog()
        expect_true(inherits(usercat, "UserCatalog"))
        expect_identical(length(usercat), 3L)
        expect_identical(urls(usercat),
            c("/api/users/user1.json",
              "/api/users/user3.json",
              "/api/users/user2.json"))
        expect_identical(names(usercat),
            c("Fake User", "Bill User", "Roger User"))
        expect_identical(emails(usercat),
            c("fake.user@example.com",
              "william.user@example.io",
              "ruser@crunch.io"))
    })
})

if (run.integration.tests) {
    with(test.authentication, {
        test_that("User can be fetched", {
            user <- try(getUser())
            expect_true(inherits(user, "ShojiObject"))
        })

        u.email <- uniqueEmail()
        u.name <- now()
        u.url <- try(invite(u.email, name=u.name, notify=FALSE))

        test_that("User can be invited", {
            skip_on_jenkins("Jenkins user needs more permissions")
            usercat <- getAccountUserCatalog()
            expect_true(u.url %in% urls(usercat))
            expect_true(u.email %in% emails(usercat))
            expect_true(u.name %in% sub(" +$", "", names(usercat)))
        })

        test_that("User can be deleted", {
            skip_on_jenkins("Jenkins user needs more permissions")
            try(crDELETE(u.url))
            usercat <- refresh(getAccountUserCatalog())
            expect_false(u.url %in% urls(usercat))
            expect_false(u.email %in% emails(usercat))
            expect_false(u.name %in% sub(" +$", "", names(usercat)))
        })

        test_that("test.user() setup/teardown", {
            skip_on_jenkins("Jenkins user needs more permissions")
            u.email <- paste0("test+", as.numeric(Sys.time()), "@crunch.io")
            u.name <- now()
            usercat <- getAccountUserCatalog()
            expect_false(u.email %in% emails(usercat))
            expect_false(u.name %in% sub(" +$", "", names(usercat)))
            with(test.user(u.email, u.name), {
                usercat <- refresh(usercat)
                expect_true(u.email %in% emails(usercat))
                expect_true(u.name %in% sub(" +$", "", names(usercat)))
                user <- index(usercat)[[u]]
                expect_false(user$account_permissions$create_datasets)
                expect_false(user$account_permissions$alter_users)
            })
            usercat <- refresh(usercat)
            expect_false(u.email %in% emails(usercat))
            expect_false(u.name %in% sub(" +$", "", names(usercat)))
        })

        test_that("User with permissions", {
            skip_on_jenkins("Jenkins user needs more permissions")
            with(test.user(advanced=TRUE), {
                user <- index(getAccountUserCatalog())[[u]]
                expect_true(user$account_permissions$create_datasets)
                expect_false(user$account_permissions$alter_users)
            })
            with(test.user(admin=TRUE), {
                user <- index(getAccountUserCatalog())[[u]]
                expect_false(user$account_permissions$create_datasets)
                expect_true(user$account_permissions$alter_users)
            })
            with(test.user(admin=TRUE, advanced=TRUE), {
                user <- index(getAccountUserCatalog())[[u]]
                expect_true(user$account_permissions$create_datasets)
                expect_true(user$account_permissions$alter_users)
            })
        })
    })

    test_that("User cannot be fetched if logged out", {
        logout()
        expect_error(getUser(),
            "You must authenticate before making this request")
    })
}
