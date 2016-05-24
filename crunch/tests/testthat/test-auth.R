context("Authentication")

test_that("On package load, the session_store exists", {
    expect_true(is.environment(session_store))
})

test_that("login checks for email and password before POSTing", {
    expect_error(crunchAuth(email=NULL),
        "Must supply the email address associated with your crunch.io account")
    expect_error(crunchAuth(email=1L, password=NULL),
        "Must supply a password")
})

with_mock_HTTP({
    test_that("Jupyter helper sets up env", {
        with(reset.option("httr_config"), {
            jupyterLogin("test_token")
            cfg <- getOption("httr_config")
            expect_identical(cfg$options$cookie, "token=test_token")
            expect_true(grepl("jupyter.crunch.io", cfg$headers[["user-agent"]]))
            expect_true(grepl("rcrunch", cfg$headers[["user-agent"]]))
        })
    })
})

if (run.integration.tests) {
    test_that("login works if crunch is running", {
        deleteSessionInfo()
        suppressMessages(login())
            expect_true("root" %in% ls(envir=session_store))
            expect_true(is.authenticated())
        logout()
        expect_false(is.authenticated())
    })

    test_that("crunchAuth succeeds when it should and not when it shouldn't", {
        logout()
        em <- getOption("crunch.email")
        pw <- getOption("crunch.pw")
        expect_true(is.character(em))
        expect_true(is.character(pw))
        expect_true(is.list(crunchAuth(em, password=pw)))
        suppressMessages(login())
        logout()
        expect_error(crunchAuth("lkjasdfksdfkjhl", password="w23nrnsod"),
            "Unable to authenticate lkjasdfksdfkjhl")
    })

    test_that("session URLs can be retrieved", {
        suppressMessages(login())
            expect_true(is.character(sessionURL("datasets")))
        logout()
        expect_error(sessionURL("datasets"),
            "You must authenticate before making this request")
    })

    test_that("login returns a session object", {
        cr <- suppressMessages(login())
            expect_true(is.list(cr))
        logout()
        expect_error(sessionURL("datasets"),
            "You must authenticate before making this request")
    })
}
