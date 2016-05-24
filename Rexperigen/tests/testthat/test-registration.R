context("registration")

## TODO: Rename context
## TODO: Add more tests
options(Rexperigen.server = "test.org")
options(Rexperigen.server.version = "2.2.0")
options(Rexperigen.experimenter = "ax")
options(Rexperigen.password = "bx")

test_that("check login", {
    expect_equal(checkLogin(), "ax")
    options(Rexperigen.experimenter = "")
    expect_error(checkLogin(), "[nN]ot logged in")
    options(Rexperigen.experimenter = "ax")
    options(Rexperigen.server.version = "1.2.0")
    expect_error(checkLogin(), "not supported")
    options(Rexperigen.server.version = "2.2.0")
})

test_that("registering new experiment", {
    with_mock(
        "RCurl::url.exists" = function(...) TRUE,
        "RCurl::getURL" = function(url, customrequest, ...){
            expect_equal(url, "http://test.org/digest/registration?experimenter=ax&sourceurl=testurl&experimentName=testname")
            expect_equal(customrequest, "POST")
            'done'
        },
        res <- registerExperiment("testurl", "testname")
    )
    expect_equal(res, "done")

})

test_that("removing registration", {
    with_mock(
        "RCurl::url.exists" = function(...) TRUE,
        "RCurl::getURL" = function(url, customrequest, ...){
            expect_equal(url, "http://test.org/digest/registration?experimenter=ax&sourceurl=testurl&experimentName=testname")
            expect_equal(customrequest, "DELETE")
            'done'
        },
        res <- removeRegistration("testurl", "testname")
    )
    expect_equal(res, "done")
})

test_that("listing registered experiments", {
    with_mock(
        "RCurl::url.exists" = function(...) TRUE,
        "RCurl::getURL" = function(url, customrequest, ...){
            expect_equal(url, "http://test.org/digest/registration?experimenter=ax")
            expect_equal(customrequest, "GET")
            '[{"sourceUrl": "testurl", "experimentName": "expname"}, {"sourceUrl": "testurl", "experimentName": "exp2"}]'
        },
        res <- getRegisteredExperiments()
    )
    expect_equal(as.character(res$sourceUrl), c("testurl","testurl"))
    expect_equal(as.character(res$experimentName), c("expname","exp2"))
})


options(Rexperigen.server = "db.phonologist.org")
options(Rexperigen.server.version = "1.0.0")
options(Rexperigen.experimenter = "")
options(Rexperigen.password = "")
