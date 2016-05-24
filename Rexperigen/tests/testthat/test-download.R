context("download")

## TODO: Rename context
## TODO: Add more tests

logoutExperigen()

options(Rexperigen.server = "test.org")
options(Rexperigen.server.version = "1.2.0")

test_that("makecsv with no auth, server 1", {
    with_mock(
        "RCurl::url.exists" = function(...) TRUE,
        "RCurl::getURL" = function(url, ...){
            expect_equal(url, "http://test.org/makecsv.cgi?sourceurl=testurl&experimentName=testname&file=default.csv")
            "column1\tcolumn2\n1\t2\n"
        },
        res <- downloadExperiment("testurl", "testname")
    )
    expect_equal(res$column1, 1)
    expect_equal(res$column2, 2)
})

options(Rexperigen.server.version = "2.2.0")

test_that("makecsv with no auth, server 2", {
    with_mock(
        "RCurl::url.exists" = function(...) TRUE,
        "RCurl::getURL" = function(url, ...){
            expect_equal(url, "http://test.org/streamresults?sourceurl=testurl&experimentName=testname&file=default.csv&ndjson=true")
            '{"column1": 1, "column2": 2}\n'
        },
        res <- downloadExperiment("testurl", "testname")
    )
    expect_equal(res$column1, 1)
    expect_equal(res$column2, 2)
})

test_that("makecsv with auth, server 2", {
    options(Rexperigen.experimenter = "aa")
    options(Rexperigen.password = "bb")
    with_mock(
        "RCurl::url.exists" = function(...) TRUE,
        "RCurl::getURL" = function(url, username, password, ...){
            expect_equal(url, "http://test.org/digest/streamresults?sourceurl=testurl&experimentName=testname&file=demographics.csv&ndjson=true")
            expect_equal(username, "aa")
            expect_equal(password, "bb")
            '{"column1": 1, "column2": 2}\n'
        },
        res <- downloadExperiment("testurl", "testname", "demographics.csv", TRUE)
    )
    expect_equal(res$column1, 1)
    expect_equal(res$column2, 2)
})


test_that("getDestinations", {
    with_mock(
        "RCurl::url.exists" = function(...) TRUE,
        "RCurl::getURL" = function(url, ...){
            expect_equal(url, "http://test.org/destinations?sourceurl=testurl&experimentName=testname")
            '["default.csv","other.csv"]'
        },
        res <- getDestinations("testurl", "testname")
    )
    expect_equal(res[1], "default.csv")
    expect_equal(res[2], "other.csv")
})

test_that("getUsers", {
    with_mock(
        "RCurl::url.exists" = function(...) TRUE,
        "RCurl::getURL" = function(url, ...){
            expect_equal(url, "http://test.org/users?sourceurl=testurl&experimentName=testname")
            'userCode\trecords\nAAA111\t23\nBBB222\t43\n'
        },
        res <- getUsers("testurl", "testname")
    )
    expect_equal(as.character(res$userCode[1]), "AAA111")
    expect_equal(res$records[2], 43)
})



options(Rexperigen.server = "db.phonologist.org")
options(Rexperigen.server.version = "1.0.0")
options(Rexperigen.experimenter = "")
options(Rexperigen.password = "")
