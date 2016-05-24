context("Teams")

with_mock_HTTP({
    test_that("Getting teams catalog", {
        teams <- try(getTeams())
        expect_true(inherits(teams, "TeamCatalog"))
        expect_identical(length(teams), 1L)
        expect_identical(names(teams), "Alpha Team")
    })

    test_that("Getting team entity", {
        teams <- try(getTeams())
        expect_true(inherits(teams[[1]], "CrunchTeam"))
        expect_true(inherits(teams$`Alpha Team`, "CrunchTeam"))
        expect_true(inherits(teams[["Alpha Team"]], "CrunchTeam"))
        expect_true(is.null(teams$`Beta Team`))
    })

    test_that("Team entity attributes", {
        ateam <- try(getTeams()[[1]])
        expect_identical(name(ateam), "Alpha Team")
        m <- try(members(ateam))
        expect_true(inherits(m, "MemberCatalog"))
        expect_identical(names(m), c("Fake User", "Roger User"))
    })
})

if (run.integration.tests) {
    with(test.authentication, {
        ucat <- getAccountUserCatalog()
        my.name <- names(ucat)[urls(ucat) == userURL()]
        my.email <- emails(ucat)[urls(ucat) == userURL()]

        teams <- try(getTeams())
        nteams.0 <- length(teams)
        test_that("Can get team catalog", {
            expect_true(inherits(teams, "TeamCatalog"))
        })

        t2 <- teams
        name.of.team1 <- now()
        test_that("Can create a team", {
            expect_false(name.of.team1 %in% names(t2))
            t2[[name.of.team1]] <- list()
            expect_true(name.of.team1 %in% names(t2))
            expect_true(length(t2) == nteams.0 + 1L)
            expect_true(inherits(t2[[name.of.team1]], "CrunchTeam"))
            expect_identical(length(members(t2[[name.of.team1]])), 1L)
            expect_identical(names(members(t2[[name.of.team1]])),
                my.name)
        })

        test_that("Can delete a team by URL", {
            t2 <- refresh(t2)
            expect_true(name.of.team1 %in% names(t2))
            try(crDELETE(self(t2[[name.of.team1]])))
            expect_false(name.of.team1 %in% names(refresh(t2)))
            ## TODO: add a delete() method for CrunchTeam, with a confirm arg.
        })

        test_that("delete method for team (requires confirmation)", {
            ## Setup
            t2 <- refresh(t2)
            nteams.2 <- length(t2)
            name.of.team2 <- now()
            expect_false(name.of.team2 %in% names(t2))
            t2[[name.of.team2]] <- list()
            expect_true(name.of.team2 %in% names(t2))
            expect_true(length(t2) == nteams.2 + 1L)

            expect_error(delete(t2[[name.of.team2]], confirm=TRUE),
                "Must confirm deleting team")
            expect_true(name.of.team2 %in% names(t2))
            expect_true(length(t2) == nteams.2 + 1L)

            ## Cleanup
            try(delete(t2[[name.of.team2]]))
            expect_false(name.of.team2 %in% names(getTeams()))
        })

        test_that("Can create a team with members", {
            skip_on_jenkins("Jenkins user needs more permissions")
            t2 <- refresh(t2)
            nteams.2 <- length(t2)
            name.of.team2 <- now()
            expect_false(name.of.team2 %in% names(t2))
            with(test.user(), {
                ucat <- getUserCatalog()
                u.email <- emails(ucat)[urls(ucat) == u]
                u.name <- names(ucat)[urls(ucat) == u]
                t2[[name.of.team2]] <- list(members=u.email)
                expect_true(name.of.team2 %in% names(t2))
                expect_true(length(t2) == nteams.2 + 1L)
                this.team <- t2[[name.of.team2]]
                expect_true(setequal(names(members(this.team)),
                    c(u.name, my.name)))
            })
            try(crDELETE(self(refresh(t2)[[name.of.team2]])))
        })

        test_that("Can add members to a team", {
            skip_on_jenkins("Jenkins user needs more permissions")
            t2 <- refresh(teams)
            name.of.team3 <- now()
            expect_false(name.of.team3 %in% names(t2))
            with(test.user(), {
                ucat <- getUserCatalog()
                u.email <- emails(ucat)[urls(ucat) == u]
                u.name <- names(ucat)[urls(ucat) == u]
                t2[[name.of.team3]] <- list()
                this.team <- t2[[name.of.team3]]
                expect_identical(names(members(this.team)),
                    my.name)
                members(this.team) <- u.email
                expect_true(setequal(names(members(this.team)),
                    c(u.name, my.name)))
            })
            try(crDELETE(self(refresh(t2)[[name.of.team3]])))
        })

        test_that("Can remove members from a team", {
            skip_on_jenkins("Jenkins user needs more permissions")
            t2 <- refresh(teams)
            name.of.team4 <- now()
            expect_false(name.of.team4 %in% names(t2))
            with(test.user(), {
                ucat <- getUserCatalog()
                u.email <- emails(ucat)[urls(ucat) == u]
                u.name <- names(ucat)[urls(ucat) == u]
                t2[[name.of.team4]] <- list()
                this.team <- t2[[name.of.team4]]
                expect_identical(names(members(this.team)),
                    my.name)
                members(this.team) <- u.email
                expect_true(setequal(names(members(this.team)),
                    c(u.name, my.name)))
                try(members(this.team)[[u.email]] <- NULL)
                expect_identical(names(members(this.team)),
                    my.name)
            })
            try(crDELETE(self(refresh(t2)[[name.of.team4]])))
        })
    })
}
