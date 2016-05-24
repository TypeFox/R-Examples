context("Projects")

with_mock_HTTP({
    projects <- session()$projects
    test_that("Getting projects catalog", {
        expect_true(inherits(projects, "ProjectCatalog"))
        expect_identical(length(projects), 1L)
        expect_identical(names(projects), "Project One")
    })

    aproject <- projects[[1]]
    test_that("Getting project from catalog", {
        expect_true(inherits(projects[[1]], "CrunchProject"))
        expect_true(inherits(projects$`Project One`, "CrunchProject"))
        expect_true(inherits(projects[["Project One"]], "CrunchProject"))
        expect_true(is.null(projects$`Beta Project`))
    })

    test_that("Project attributes", {
        expect_identical(name(aproject), "Project One")
    })

    test_that("Simple project creation by assignment", {
        expect_error(projects[["A new project"]] <- list(),
            'POST /api/projects.json {"name":"A new project"}',
            fixed=TRUE)
    })

    test_that("Project deletion", {
        expect_error(delete(projects[[1]], confirm=TRUE),
            "Must confirm deleting project")
        with(consent(), expect_error(delete(projects[[1]], confirm=TRUE),
            "DELETE /api/projects/project1.json"))
    })

    m <- members(aproject)
    test_that("Project members catalog", {
        expect_true(inherits(m, "MemberCatalog"))
        expect_identical(names(m), c("Fake User", "Roger User"))
    })

    test_that("Add members by members<-", {
        expect_error(members(aproject) <- c("new.user@crunch.io", "foo@example.co"),
            'PATCH /api/projects/project1/members.json {"new.user@crunch.io":{},"foo@example.co":{}}',
            fixed=TRUE)
    })

    test_that("Add members doesn't re-add if already a member", {
        expect_error(members(aproject) <- c("new.user@crunch.io", "roger.user@example.com"),
            'PATCH /api/projects/project1/members.json {"new.user@crunch.io":{}}',
            fixed=TRUE)
    })

    test_that("Remove members by <- NULL", {
        expect_error(members(aproject)[["roger.user@example.com"]] <- NULL,
            'PATCH /api/projects/project1/members.json {"roger.user@example.com":null}',
            fixed=TRUE)
    })

    d <- datasets(aproject)
    test_that("Project datasets catalog", {
        expect_true(inherits(d, "DatasetCatalog"))
        expect_identical(names(d), "ECON.sav")
    })

    test_that("Add datasets to project by <- a dataset (which transfers ownership)", {
        ds <- loadDataset("test ds")
        expect_error(datasets(aproject) <- ds,
            'PATCH /api/datasets/dataset1.json {"owner":"/api/projects/project1.json"}',
            fixed=TRUE)
    })
})

if (run.integration.tests) {
    with(test.authentication, {
        projects <- session()$projects
        ucat <- getAccountUserCatalog()
        my.name <- names(ucat)[urls(ucat) == userURL()]
        my.email <- emails(ucat)[urls(ucat) == userURL()]

        nprojects.0 <- length(projects)
        test_that("Can get project catalog", {
            expect_true(inherits(projects, "ProjectCatalog"))
        })

        t2 <- projects
        name.of.project1 <- now()
        test_that("Can create a project", {
            expect_false(name.of.project1 %in% names(t2))
            t2[[name.of.project1]] <- list()
            expect_true(name.of.project1 %in% names(t2))
            expect_true(length(t2) == nprojects.0 + 1L)
            expect_true(inherits(t2[[name.of.project1]], "CrunchProject"))
            expect_identical(length(members(t2[[name.of.project1]])), 1L)
            expect_identical(names(members(t2[[name.of.project1]])),
                my.name)
        })

        test_that("Can delete a project by URL", {
            t2 <- refresh(t2)
            expect_true(name.of.project1 %in% names(t2))
            try(crDELETE(self(t2[[name.of.project1]])))
            expect_false(name.of.project1 %in% names(refresh(t2)))
        })

        test_that("Can create a project with members", {
            skip("TODO")
            skip_on_jenkins("Jenkins user needs more permissions")
            t2 <- refresh(t2)
            nprojects.2 <- length(t2)
            name.of.project2 <- now()
            expect_false(name.of.project2 %in% names(t2))
            with(test.user(), {
                ucat <- getUserCatalog()
                u.email <- emails(ucat)[urls(ucat) == u]
                u.name <- names(ucat)[urls(ucat) == u]
                t2[[name.of.project2]] <- list(members=u.email)
                expect_true(name.of.project2 %in% names(t2))
                expect_true(length(t2) == nprojects.2 + 1L)
                this.project <- t2[[name.of.project2]]
                expect_true(setequal(names(members(this.project)),
                    c(u.name, my.name)))
            })
            try(crDELETE(self(refresh(t2)[[name.of.project2]])))
        })

        test_that("Can add members to a project", {
            skip_on_jenkins("Jenkins user needs more permissions")
            with(cleanup(testProject()), as="tp", {
                with(test.user(), {
                    ucat <- getUserCatalog()
                    u.email <- emails(ucat)[urls(ucat) == u]
                    u.name <- names(ucat)[urls(ucat) == u]
                    expect_identical(names(members(tp)),
                        my.name)
                    members(tp) <- u.email
                    expect_true(setequal(names(members(tp)),
                        c(u.name, my.name)))
                })
            })
        })

        test_that("Can remove members from a project", {
            skip_on_jenkins("Jenkins user needs more permissions")
            with(cleanup(testProject()), as="tp", {
                with(test.user(), {
                    ucat <- getUserCatalog()
                    u.email <- emails(ucat)[urls(ucat) == u]
                    u.name <- names(ucat)[urls(ucat) == u]
                    expect_identical(names(members(tp)),
                        my.name)
                    members(tp) <- u.email
                    expect_true(setequal(names(members(tp)),
                        c(u.name, my.name)))
                    try(members(tp)[[u.email]] <- NULL)
                    expect_identical(names(members(tp)),
                        my.name)
                })
            })
        })

        test_that("Can add datasets to projects", {
            with(test.dataset(df), {
                with(cleanup(testProject()), as="tp", {
                    expect_true(inherits(tp, "CrunchProject"))
                    expect_identical(length(datasets(tp)), 0L)
                    datasets(tp) <- ds
                    expect_identical(names(datasets(tp)), name(ds))
                })
            })
        })
    })
}
