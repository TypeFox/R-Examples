context("Server functions")

test_that("setExperigenServer", {
    with_mock(
        "RCurl::url.exists" = function(...) {FALSE},
        {
            expect_error(setExperigenServer("test.server"),
                         "No server")
        }
    )
    with_mock(
        "RCurl::getURL" = function(url, ...){
            expect_equal(url, "http://test.server/version")
            "3.1.2"
        },
        "RCurl::url.exists" = function(...) {TRUE},
        {
            setExperigenServer("test.server")
            expect_equal(getOption("Rexperigen.server"), "test.server")
            expect_equal(getOption("Rexperigen.server.version"), "3.1.2")
        },
        setExperigenServer
    )
    expect_equal(versionMain(), 3)
})

test_that("Setting credentials", {
    setExperigenCredentials("a", "bb", check = FALSE, quiet = TRUE)
    expect_equal(getOption("Rexperigen.experimenter"), "a")
    expect_equal(getOption("Rexperigen.password"), "bb")
    
    ## server is test.server now
    with_mock(
        "RCurl::url.exists" = function(...) TRUE,
        "RCurl::getURL" = function(url, username, password, ...){
            expect_equal(url, "http://test.server/digest/me")
            expect_equal(username, "apple")
            expect_equal(password, "pear")
            "apple"
        },
        {
            val <- setExperigenCredentials("apple", "pear", quiet = TRUE)
            expect_equal(getOption("Rexperigen.experimenter"), "apple")
            expect_equal(getOption("Rexperigen.password"), "pear")
            expect_true(val)
        }
    )
    with_mock(
        "RCurl::url.exists" = function(...) TRUE,
        "RCurl::getURL" = function(url, username, password, ...){
            expect_equal(url, "http://test.server/digest/me")
            expect_equal(username, "apple")
            expect_equal(password, "pear")
            "Unauthorized"
        },
        {
            val <- setExperigenCredentials("apple", "pear", quiet = TRUE)
            expect_equal(getOption("Rexperigen.experimenter"), "")
            expect_equal(getOption("Rexperigen.password"), "")
            expect_false(val)
        }
    )
    options(Rexperigen.server.version = "1.0.0")
    expect_error({
        val <- setExperigenCredentials("x", "y")
    }, "not supported")
    expect_false(val)
    expect_equal(getOption("Rexperigen.experimenter"), "")
    expect_equal(getOption("Rexperigen.password"), "")
    options(Rexperigen.server.version = "3.1.2")
})

test_that("creating experimenters", {
    with_mock(
        "RCurl::url.exists" = function(...) TRUE,
        "RCurl::getURL" = function(url, customrequest, ...){
            if(customrequest == "POST"){
                ha1 <- digest::digest("abc:Experimenters:def", algo="md5", serialize = FALSE)
                expect_equal(url, paste0("http://test.server/experimenter?experimenter=abc&ha1=", ha1))
                return("done")
            }
            else {
                expect_equal(url, paste0("http://test.server/digest/me"))
                return("abc")
            }
        },
        createExperimenter("abc", "def")
    )
    expect_equal(getOption("Rexperigen.experimenter"), "abc")
    expect_equal(getOption("Rexperigen.password"), "def")
    options(Rexperigen.server.version = "1.0.0")
    expect_error({
        createExperimenter("x", "y")
    }, "not supported")
    options(Rexperigen.server.version = "3.1.2")
})


options(Rexperigen.server = "db.phonologist.org")
options(Rexperigen.server.version = "1.0.0")
